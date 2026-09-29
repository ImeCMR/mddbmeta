#!/bin/bash
set -e
G=/usr/local/gromacs/avx2_256/bin/gmx
cd "$(dirname "$0")/prep"
export GMX_MAXBACKUP=-1
$G editconf -f ala_ala.gro -o box.gro -bt dodecahedron -d 1.0 -c >>../build.log 2>&1
$G solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro >>../build.log 2>&1
cat > em.mdp <<M
integrator = steep
emtol = 1000
nsteps = 500
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 0.9
rvdw = 0.9
pbc = xyz
M
cat > nvt.mdp <<M
integrator = md
dt = 0.002
nsteps = 1000
nstxout-compressed = 250
nstenergy = 250
nstlog = 250
continuation = no
gen_vel = yes
gen_temp = 300
gen_seed = 11
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 0.9
rvdw = 0.9
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = 300 300
pcoupl = no
pbc = xyz
M
sed -e 's/nsteps = 1000/nsteps = 1000/' -e 's/continuation = no/continuation = yes/' -e 's/gen_vel = yes/gen_vel = no/' -e 's/pcoupl = no/pcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\ncompressibility = 4.5e-5/' nvt.mdp > npt.mdp
sed -e 's/nsteps = 1000/nsteps = 2500/' -e 's/nstxout-compressed = 250/nstxout-compressed = 500/' npt.mdp > ../md.mdp
$G grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 5 >>../build.log 2>&1
$G mdrun -deffnm em -nt 4 -nb cpu >>../build.log 2>&1
$G grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 5 >>../build.log 2>&1
$G mdrun -deffnm nvt -nt 4 -nb cpu >>../build.log 2>&1
$G grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 5 >>../build.log 2>&1
$G mdrun -deffnm npt -nt 4 -nb cpu >>../build.log 2>&1
# rep1: three chained segments, each a new tpr continuing from the previous checkpoint
cd ../rep1
cp ../md.mdp md.mdp
$G grompp -f md.mdp -c ../prep/npt.gro -t ../prep/npt.cpt -p ../prep/topol.top -o md_1.tpr -maxwarn 5 >>../build.log 2>&1
$G mdrun -deffnm md_1 -nt 4 -nb cpu >>../build.log 2>&1
for i in 2 3; do p=$((i-1))
  $G grompp -f md.mdp -c md_$p.gro -t md_$p.cpt -p ../prep/topol.top -o md_$i.tpr -maxwarn 5 >>../build.log 2>&1
  $G mdrun -deffnm md_$i -nt 4 -nb cpu >>../build.log 2>&1
done
# rep2: one tpr, extended twice via checkpoint continuation (append, then -noappend)
cd ../rep2
cp ../md.mdp md.mdp
$G grompp -f md.mdp -c ../prep/npt.gro -t ../prep/npt.cpt -p ../prep/topol.top -o prod.tpr -maxwarn 5 >>../build.log 2>&1
$G mdrun -deffnm prod -nt 4 -nb cpu >>../build.log 2>&1
$G convert-tpr -s prod.tpr -extend 5 -o prod.tpr >>../build.log 2>&1
$G mdrun -deffnm prod -cpi prod.cpt -nt 4 -nb cpu >>../build.log 2>&1
$G convert-tpr -s prod.tpr -extend 5 -o prod.tpr >>../build.log 2>&1
$G mdrun -deffnm prod -cpi prod.cpt -noappend -nt 4 -nb cpu >>../build.log 2>&1
echo BUILD_OK
