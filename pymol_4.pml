
# Example: investigate possible interaction
# involving Lys11 in ubiquitin (1UBQ)



reinitialize


fetch 1ubq, async=0


remove solvent


hide everything
show cartoon, 1ubq
bg_color white



show sticks, 1ubq and resi 11
zoom 1ubq and resi 11, 8

# Label residue identity
label 1ubq and resi 11 and name CA, "%s%s" % (resn,resi)

# Label all atoms
label 1ubq and resi 11, name

# For Lys we should see:
# N, CA, C, O, CB, CG, CD, CE, NZ
# NZ is the terminal side-chain nitrogen



# Find residues close to Lys11


# Find complete residues within 4 Å of Lys11
select nearby, byres (1ubq within 4 of (1ubq and resi 11))

show sticks, nearby

# Label nearby residues by name and number
label nearby and name CA, "%s%s" % (resn,resi)

# Now visually inspect for possible interacting residues,
# especially ASP or GLU for a possible salt bridge.


# Suppose Glu34 is identified as a candidate


show sticks, 1ubq and resi 34

# Label its atoms
label 1ubq and resi 34, name

# Glu should contain:
# CB, CG, CD, OE1, OE2
# OE1/OE2 are the carboxylate oxygens.



# Ask PyMOL to identify possible polar contacts


distance hbonds, 1ubq and resi 11, 1ubq and resi 34, mode=2

# mode=2 searches for polar / H-bond-like contacts



# Measure the exact relevant atom distances


distance d1, 1ubq and resi 11 and name NZ, \
             1ubq and resi 34 and name OE1

distance d2, 1ubq and resi 11 and name NZ, \
             1ubq and resi 34 and name OE2


# Print exact numerical distances in the PyMOL console
get_distance 1ubq and resi 11 and name NZ, \
             1ubq and resi 34 and name OE1

get_distance 1ubq and resi 11 and name NZ, \
             1ubq and resi 34 and name OE2



# Make the result easier to see


set dash_color, black
zoom 1ubq and resi 11+34, 6

# Optional: clean labels and only label residues
hide labels
label 1ubq and resi 11+34 and name CA, "%s%s" % (resn,resi)