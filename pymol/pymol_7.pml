reinitialize
fetch 1stp, async=0
remove solvent
hide everything
show cartoon, 1stp
color grey80, 1stp
select ligand, 1stp and resn BTN
# stp (streptavidin) is the proteins
# btn (biotin) is the ligand
show sticks, ligand
util.cbac ligand

select pocket, (byres (ligand around 4)) and polymer
# pocket for binding pocket
show sticks, pocket and not name N+C+O
get_area, pocket
util.cbay pocket
label pocket and name CA, "%s%s" % (resn, resi)

distance hbonds, ligand, pocket, mode=2
# finding polar contacts between ligand and pocket
set dash_color, black
iterate pocket and name CA, print(resn, resi)
zoom ligand, 5

