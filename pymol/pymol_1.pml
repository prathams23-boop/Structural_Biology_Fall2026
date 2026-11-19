reinitialize

fetch 1ubq, async=0
remove solvent
remove inorganic
hide everything
dss 
# assigned secondary structure
show cartoon, 1ubq

set_cartoon_fancy_helices, on
color grey80, 1ubq
color skyblue, 1ubq and ss H
color salmon, 1ubq and ss S 

bg_color white
orient 1ubq
count_atoms 1ubq and name CA
