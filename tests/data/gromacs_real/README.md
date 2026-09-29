# Real GROMACS fixture

Generated with GROMACS 2022.3 (NGC container, CPU) by `build.sh`, starting from
MDDB-workflow's dummy Ala-Ala files (`test/data/input/dummy/gromacs`, Apache-2.0):
solvated in a rhombic dodecahedron (TIP3P, amber99sb-ildn), then

- `prep/`: em (steep) -> nvt -> npt (`grompp -t` continuations);
- `rep1/`: three production segments chained with `grompp -t md_N.cpt` (each restarts the clock);
- `rep2/`: one production run extended twice with `mdrun -cpi` (appended), then `-noappend`.

Trimmed to the files mddbmeta reads (no .cpt/.edr/.trr, one .gro). Working directories and
host names in the logs were replaced with placeholders.
