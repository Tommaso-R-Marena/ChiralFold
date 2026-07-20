# ChiralFold mirror-image comparison template.
# Replace l_structure.pdb and d_structure.pdb with your structures.

reinitialize
bg_color white
load l_structure.pdb, l_structure
load d_structure.pdb, d_structure
hide everything, l_structure or d_structure
show cartoon, l_structure or d_structure
show sticks, (l_structure or d_structure) and not elem H
color marine, l_structure
color magenta, d_structure
translate [35, 0, 0], d_structure
label l_structure and name CA and resi 1, "L structure"
label d_structure and name CA and resi 1, "D/mirror structure"
set cartoon_transparency, 0.08
set ray_opaque_background, off
group mirror_comparison, l_structure d_structure
zoom l_structure or d_structure
