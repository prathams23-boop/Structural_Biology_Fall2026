# distance bn lys11 and glu34 in ubq

reinitialize
fetch 1ubq, async=0
remove solvent
hide everything
show cartoon, 1ubq
show sticks, 1ubq and resi 11+34 and not name N+C+O 

# command for distance
# distance name, atom1, atom2

distance d1, 1ubq and resi 11 and name NZ, 1ubq and resi 34 and name OE1
distance d2, 1ubq and resi 11 and name NZ, 1ubq and resi 34 and name OE2
get_distance 1ubq and resi 11 and name NZ, 1ubq and resi 34 and name OE2
get_distance 1ubq and resi 11 and name NZ, 1ubq and resi 34 and name OE1

set dash_color, black
label 1ubq and resi 11+34 and name CA, "%s%s" % (resn,resi)
zoom 1ubq and resi 11+34, 6

# if we forget residue structure
 label 1ubq and resi 34, name
 label 1ubq and resi 11, name
 bg_color white