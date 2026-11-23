# Changelog

## Unreleased: MELD topology

- `mddbmeta meld-prmtop <run>` writes an AMBER prmtop (`topology.prmtop`) from a MELD run's `Data/system.dat`:
  - read with an inert unpickler that imports nothing (no MELD needed, no code from the file runs);
  - Grappa's split bonded forces are merged and its zeroed ff14SB torsions dropped, so the prmtop holds the simulated parameters;
  - implicit solvent is left out (no prmtop form); any other force a prmtop cannot hold is refused;
  - the prmtop's bond/angle/torsion/nonbonded energies must match the MELD system at the saved coordinates, or nothing is written.
- `openmm` and `parmed` (and `numpy`) are now dependencies.
- `MINE004` on a MELD run without a topology points to `mddbmeta meld-prmtop`.
- Discovery writes `topology.prmtop` from each MELD run's `Data/system.dat` before scanning, and MELD runs use it as their topology (other prmtops are listed as not used). Reused while newer than `system.dat`; never overwrites a file it did not write; `--no-meld-prmtop` keeps discovery read-only. New findings: `MELD005` (topology written/reused), `MELD006` (could not be written).
- MELD setup scripts in the older syntax (`SystemBuilder(forcefield=...)`, `options = RunOptions(); options.use_bigger_timestep = True`, ...) are read: force field, implicit-solvent model, timestep, steps per exchange and minimisation steps. Before, these were missed and the timestep fell back to 2 fs. `AmberOptions(forcefield=...)` now gives `ff` too, `ff14sbside` is reported as `Amber ff14SBonlysc` (the leaprc MELD sources), and when both `use_big_timestep` and `use_bigger_timestep` are set the timestep is 3.5 fs, as in MELD.
- MELD <= 0.4's single `remd.log` (all ranks interleaved behind a `host:rank date time LEVEL logger:` prefix) is read: program, MELD and OpenMM_Meld versions, restarts and progress. Before, only `Logs/remd_NNN.log` was recognised, so program and version were left unstated.
- Much faster discovery of long MELD runs (1kh0: 2470 blocks, 246900 frames, from about 50 minutes to a few):
  - the data store is read one block's own frames at a time. MELD's blocks share one growing time axis, so reading each block's whole `alphas` variable cost more with every block, and the per-value Python loop is now numpy;
  - a PDB is read up to its first `END`/`ENDMDL`: MELD's multi-frame `Data/trajectory.pdb` (13.9 GB on 1kh0) was read to the end;
  - in a numbered series of files (`block_000001.nc`, ...), once three members were identified from their own bytes, the rest are identified by name without opening them.
- Fewer, clearer warnings: runs that share one log and did not finish (every replica of a REMD job) give one `CMP002` with the log's reason (e.g. "stopped at exchange 246992 of 640000") instead of one per run; an unreadable MELD data store gives one `LOAD001` saying temperatures are left blank, instead of one `MELD001` per replica that wrongly said the initial ladder was used.
- `MELD009` (warning): ladder DCDs of different lengths, i.e. `extract_traj_dcd` runs that stopped early.
- A starting-structure PDB that mddbmeta wrote under an older name (`meld_system.pdb`) is never taken for the run's own PDB.
- MELD runs with no PDB of their own get their starting structure too: `topology.pdb` is written from the template coordinates in `Data/system.dat` (same atoms and order as the prmtop) and used as `input_structure_filepath` when no PDB in the run matches the topology. `meld-prmtop` writes it as well (`--no-pdb` to skip). Older MELD `system.dat` files (`meld.system.system.System`, which keep the AMBER prmtop and inpcrd as text) are supported too: the saved prmtop is written out unchanged apart from its title, and the starting structure comes from the saved inpcrd (angstrom). New findings: `MELD007` (structure written/reused and used), `MELD008` (could not be written).
- Adapters get a `prepare(root, max_depth, exclude)` hook that runs before the scan.

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
