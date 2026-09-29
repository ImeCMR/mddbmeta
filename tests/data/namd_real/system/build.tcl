package require psfgen
topology ../toppar/top_all36_prot.rtf
topology ../toppar/toppar_water_ions_namd.str
pdbalias atom ALA OXT OT2
segment PROA { pdb ala_heavy.pdb }
coordpdb ala_heavy.pdb PROA
guesscoord
writepsf pep.psf
writepdb pep.pdb
package require solvate
solvate pep.psf pep.pdb -t 8 -o solv
mol delete all
mol new solv.psf; mol addfile solv.pdb
set all [atomselect top all]
set mm [measure minmax $all]
set c [measure center $all]
set d [vecsub [lindex $mm 1] [lindex $mm 0]]
set fh [open cell.txt w]
puts $fh "[lindex $d 0] [lindex $d 1] [lindex $d 2] $c"
close $fh
# restraint reference: heavy protein atoms get B = 1
$all set beta 0
[atomselect top "segname PROA and noh"] set beta 1
$all writepdb restraint.pdb
puts "NATOMS [$all num]"
exit
