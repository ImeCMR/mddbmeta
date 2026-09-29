# Real OpenMM fixture

Generated with OpenMM 8.6.1 (CPU) from MDDB-workflow's dummy Ala-Ala PDB (Apache-2.0), solvated
with `Modeller.addSolvent` (amber14-all + amber14/tip3pfb, 0.15 M NaCl):

- `prep/setup.py`: minimize, NVT then NPT in one script (the barostat is added between them);
  serializes `system.xml`, `integrator.xml` and `eq_state.xml`.
- `rep1/run.py`: one script looping over three segments with f-string file names; each segment
  after the first loads `prod_{i-1}.chk`.
- `rep2/prod.py` (XTC reporter, new velocities from `eq_state.xml`) and `rep2/restart.py`
  (continues from `prod.chk`).

Binary checkpoints (`*.chk`) are omitted: mddbmeta only needs their names from the scripts.
