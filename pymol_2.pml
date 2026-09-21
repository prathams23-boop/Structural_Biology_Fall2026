reinitialize

fetch 1ubq, async=0
remove solvent
remove inorganic
hide everything
show cartoon, 1ubq
color white, 1ubq

# selecting residues

select patch, 1ubq and resi 8+44+70
# this is ubqs canonical hydrophobic patch
# This surface patch serves as the primary recognition site for proteasomes, ubiquitin-binding domains (UBDs), and shuttling proteins
show sticks, patch and not name N+C+O 
color orange, patch

select acidic, 1ubq and resn ASP+GLU
select basic, 1ubq and resn LYS+ARG
color red, acidic
color blue, basic

select near44, (byres((1ubq and resi 44) around 5)) and polymer
show lines, near44
label patch and name CA, "%s%s" % (resn,resi)

set label_size, 18
iterate near44 and name CA, print(resn,resi)
count_atoms basic and name CA
count_atoms acidic and name CA
zoom patch, 5

