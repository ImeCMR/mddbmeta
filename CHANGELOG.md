# Changelog

## Unreleased: M6 MELD

- MELD adapter:
  - setup script (static): run constants, RunOptions, Grappa/Amber builder options → timestep rule, temperature scaler, adaptor;
  - leader/worker `remd` logs: versions, sessions/restarts, completion;
  - data-store `alphas` (pickles never opened);
  - `extract_trajectory` scripts → ladder vs walker DCDs.
- Ladder positions are exported as `type: ensemble` MDs, with per-MD temperature from the final alpha and `metadditions.meld`.
- New findings: MELD001 (adapted ladder), MELD002 (extracted vs stored frames), MELD003 (walkers not exported), MELD004 (restarts).
- Core:
  - adapters can declare replicas (`declare_replicas`);
  - adapter-set roles are kept;
  - `type: ensemble` suppresses `framestep`;
  - `method`/`metadditions` are exported when mined;
  - replicas sharing a directory get their id as `mdir`.

## Unreleased: M5 OpenMM

- OpenMM adapter:
  - static analysis of run scripts (`ast`, never executed): units, integrators, thermostats, barostats, `createSystem` options, reporters, `step()`, load/save, `setTemperature` ramps, restraints, `Modeller`, f-string file templates with placeholder arithmetic;
  - StateDataReporter logs, including stdout mixed with other output;
  - `State`/`System`/`Integrator` XML.
- Steps are StateDataReporter segments, matched to scripts by file template. Continuations come from load/save pairs or from the same script continuing.
- OpenMM DCD/XTC frame times are reporter-relative and are now shifted onto the run clock.
- The log timestep overrides the script's.
- Structure-only systems export `input_topology_filepath: 'no'` (MDDB-workflow's no-topology flag, `MINE004`).
- New core hooks: `run_evidence`, `BuildEvidence.linked`; PDB residue parsing; PDB never counts as another engine for `ENG001`.
- Real fixtures: `tests/data/openmm_real`, `tests/data/openmm_amber`. `tests/data` is excluded from ruff.

## Unreleased: M4 NAMD

- NAMD adapter:
  - config reader (`set`/`$var`/`source`; Tcl blocks recorded, not evaluated);
  - log reader (`Info:` resolved settings, `TCL:` run/minimize counts and parameter changes, `ENERGY:` rows, FATAL detection);
  - dcd, binary `.coor`/`.vel`, xsc, psf and pdb readers.
- Exact continuations from `BINARY COORDINATES` = another run's `.coor` or `.restart.coor`.
- Heating ramps from logged parameter changes.
- CHARMM force-field labels, gated by the residue classes in the psf.
- `amber on` runs:
  - the adapter adopts the AMBER prmtop/inpcrd (new `claim_foreign` hook);
  - `ENG001` now only counts files the chosen engine does not use.
- `input_structure_filepath` exported when all production runs read one coordinate file (NAMD's `COORDINATE PDB`).
- Adapter-resolved continuations no longer make the starting structure look ambiguous.
- Real fixtures: `tests/data/namd_real` (NAMD 2.14b1, CHARMM36) and `tests/data/namd_amber`.
- `MRG001` corrected: MDDB-workflow merges parts with MDTraj `mdconvert` (list order, nothing dropped), not `gmx trjcat`. The check now reports frames duplicated at joins, as a warning.

## Unreleased: M3 GROMACS

- GROMACS adapter, read with the stdlib only (no GROMACS needed):
  - mdp;
  - mdrun logs: multiple sessions, checkpoint reads, resolved Input Parameters, energy records;
  - gro and top;
  - xtc/trr frame walking; tpr header (version, atom count).
- Continuations from `mdrun -cpi` checkpoint reads and from energy fingerprints (`LIN003`).
- New checks:
  - `CONT005`: clock restarted on a continuation;
  - `MRG001`: frames duplicated at joins by MDDB-workflow's merge (a warning). An earlier draft assumed mwf merges xtc with `gmx trjcat` and would drop frames of clock-reset parts. mwf actually selects MDTraj `mdconvert` first, which keeps every frame in list order (verified on real files), so the check was corrected before release.
- The core now knows GROMACS writes the step-0 frame (frame counts, trajectory joins, start times).
- Box shapes: Dodecahedron and GROMACS's triclinic truncated octahedron.
- Replica directories whose run names don't overlap are still replicas when no restart links them. Linked chunk directories stay one replica.
- `ENG001` (mixed engines) is now a warning; the engine with the most run logs wins.
- Real GROMACS 2022.3 fixture in `tests/data/gromacs_real`.

## 0.1.0 (unreleased): M0 + M1

- Engine-neutral model (Project → Replica → Phase → Step) with per-value provenance.
- AMBER adapter:
  - prmtop, mdin namelists and mdout (File Assignments, CONTROL DATA, completion);
  - ASCII rst7/mdcrd and NetCDF trajectories/restarts;
  - tleap scripts and logs.
- Discovery:
  - content sniffing;
  - replica inference, refusing ambiguous layouts;
  - one role classifier;
  - restart-lineage chains and production-chain ordering.
- Validation:
  - continuity (restart and trajectory joins);
  - sequence holes;
  - completeness;
  - consistency.
- Mining of MDDB fields:
  - program, version, timestep, framestep, temp, ensemble, boxtype, wat, ff, type;
  - `mds[]` and topology.
- Export:
  - MDDB inputs plus a provenance sidecar;
  - `--fill`, `reconcile`, and a `project.yaml` manifest with overrides;
  - generator and template for `mwf dataset inputs`.
- CLI: `info`, `discover`, `validate`, `export`, `reconcile`, `template-path`, `generator-path`.
