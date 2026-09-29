# Real NAMD fixture

Generated with NAMD 2.14b1 (verbs-smp, CPU) on a solvated tri-alanine built with VMD 1.9.4 psfgen +
solvate (CHARMM36 protein, CHARMM TIP3P), starting from MDDB-workflow's dummy Ala-Ala PDB (Apache-2.0):

- `prep/`: minimize -> heating (Tcl loop raising langevinTemp 50 -> 300 K) -> NPT with harmonic restraints;
- `rep1/`: three segments, each reading the previous run's final `.coor/.vel/.xsc` with `firsttimestep`;
- `rep2/`: `prod` (Tcl variables) and `prod_2`, which restarts from `prod.restart.*` and reads its
  first timestep from the `.xsc` with a Tcl proc.

Trimmed to what mddbmeta reads; `heat.dcd` and `eq.dcd` are omitted, so discovery reports them missing (CMP001). The CHARMM parameter files (from VMD's bundle) are not included;
the logs still name them. Working directories, host and user names in the logs were replaced.
`system/build.tcl` is the psfgen/solvate script.
