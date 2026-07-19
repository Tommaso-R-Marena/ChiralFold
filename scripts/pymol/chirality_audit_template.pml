# ChiralFold chirality audit template.
# Replace input.pdb and edit chirality_wrong to the residues reported by ChiralFold.

reinitialize
bg_color white
load input.pdb, audited
hide everything, audited
show cartoon, audited
show sticks, audited and not elem H
color forest, audited
color gray70, audited and resn GLY

# Example: select chirality_wrong, (audited and chain A and resi 3) or (audited and chain A and resi 8)
select chirality_wrong, none
color red, chirality_wrong
show spheres, chirality_wrong and name CA+CB
label chirality_wrong and name CA, "%s%s" % (resn, resi)

set sphere_scale, 0.35
set ray_opaque_background, off
zoom audited
