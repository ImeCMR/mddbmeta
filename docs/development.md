# Developing mddbmeta

Notes for contributors. Users only need the [README](../README.md).

## Architecture

```
src/mddbmeta/
  provenance.py   Mined(value, unit, sources, method, confidence), Finding, FileLoadError
  findings.py     stable finding codes (catalogue)
  model.py        Project / Replica / Phase / Step, StepSettings / StepRuntime / TrajectoryInfo
  engines/
    base.py       EngineAdapter: sniff, parse, build_steps, system_info, coords_natom_time, ...
    registry.py   built-ins + the `mddbmeta.engines` entry-point group
    amber/        prmtop, mdin (namelists), mdout, rst7/mdcrd, NetCDF, tleap logs, semantics
    gromacs/      mdp, mdrun log (sessions), gro/top, xtc/trr/tpr/cpt (stdlib XDR), semantics
    namd/         config (set/$var/source), log (Info/TCL/ENERGY), dcd/coor/xsc/psf/pdb, semantics
    openmm/       run-script static analysis (ast), StateDataReporter logs, XmlSerializer files
    meld/         setup-script analysis, remd logs, data-store alphas, extract_trajectory scripts
  discovery/      scan (content sniffing), roles, replicas, lineage (input coords → chains)
  validate/       continuity, completeness, consistency, schema_lite
  mining/         fields (MDDB rules), water, boxtype, forcefield
  export/         mddb (inputs + sidecar), reconcile (fill/compare), manifest (project.yaml)
  integrations/   mwf_generator.py + Jinja template for `mwf dataset inputs`
  cli/            argparse front end and text rendering
  pipeline.py     discover(): scan → steps → roles → replicas → lineage → validate → mine
```

Design rules:

- Discovery, validation, mining and export are **engine-neutral**. Only `engines/<name>/` knows file formats.
- Parsers never raise: failures become `FileLoadError`s and `LOAD001` findings.
- Every inference is announced as a finding.
- Paths on the model are POSIX, relative to the project root. The root is resolved to an absolute path on entry.

## Adding an engine

1. Implement `EngineAdapter` in `engines/<engine>/`. `sniff` classifies files by content. `parse` never raises. `build_steps` groups files into `Step`s and fills `settings`, `runtime`, `trajectory` and the `FileRole.INPUT_COORDS` / `OUTPUT_RESTART` files; lineage and everything downstream is derived from those.
2. Register it under the `mddbmeta.engines` entry point, or in `engines/registry.py`.
3. Add a synthetic-file writer and scenarios under `tests/`.

## Tests and lint

```bash
module load python/3.11                      # on HiPerGator
python -m venv .venv && .venv/bin/pip install -e ".[dev,yaml,netcdf]"
.venv/bin/pytest                             # synthetic fixtures; no third-party data vendored
.venv/bin/pytest --update-golden             # after an intended output change
MDDBMETA_EXTERNAL_DATA=/path/to/runs .venv/bin/pytest -m external_data
.venv/bin/ruff check src tests && .venv/bin/ruff format --check src tests
```

The AMBER tests generate valid synthetic files (`tests/amberfiles.py`).

The GROMACS tests use `tests/data/gromacs_real/`: genuine GROMACS 2022.3 output, reproducible with its `build.sh`.

The OpenMM tests use `tests/data/openmm_real/` (OpenMM 8.6.1: a Modeller-built system, a looped segment script and a checkpoint restart) and `tests/data/openmm_amber/` (a prmtop run logging to stdout).

Fixture inputs under `tests/data/` are excluded from ruff so they stay byte-identical to what ran.

The NAMD tests use `tests/data/namd_real/` (NAMD 2.14b1, CHARMM36, built with VMD psfgen) and `tests/data/namd_amber/` (an `amber on` run). Their READMEs describe how each was made. They also read MDDB-workflow's own `test/data/input` files when present.

After changing finding codes, run `python scripts/gen_findings_doc.py`.

