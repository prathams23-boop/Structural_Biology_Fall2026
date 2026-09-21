# bond angle
# the N-Calpha-C angle of Gly10

reinitialize

fetch 1ubq, async=0
remove solvent
hide everything

show sticks, resi 9-11
label 1ubq and resi 9, name
label 1ubq and resi 10, name
label 1ubq and resi 11, name
label 1ubq and resi 9 and name CA, "%s%s" % (resn,resi)
label 1ubq and resi 10 and name CA, "%s%s" % (resn,resi)
label 1ubq and resi 11 and name CA, "%s%s" % (resn,resi)

angle a1, resi 10 and name N, resi 10 and name CA, resi 10 and name C
get_angle, resi 10 and name N, resi 10 and name CA, resi 10 and name C

zoom resi 10,4

