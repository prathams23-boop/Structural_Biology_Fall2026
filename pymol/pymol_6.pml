# measuring torsion angles, for gly10 and phe45

reinitialize
fetch 1ubq, async=0
remove solvent 
hide everything
show sticks, resi 9-11
# phi of resi 10=C(9)-N(10)-CA(10)-C(10)

dihedral phi10, resi 9 and name C, resi 10 and name N, resi 10 and name CA, resi 10 and name C

# omega between 10 and 11 = CA(10) - C(10) - N(11) - CA(11)

get_dihedral resi 10 and name CA, resi 10 and name C, resi 11 and name N, resi 11 and name CA

# chi1 of residue 45 = N - CA - CB - CG

get_dihedral resi 45 and name N, resi 45 and name CA, resi 45 and name CB, resi 45 and name CG

#shortcut

phi_psi resi 10


