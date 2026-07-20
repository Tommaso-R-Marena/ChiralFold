# ChiralFold AF3 correction before/after template.
# Replace before.pdb and after.pdb, then edit corrected_before/after selections.

reinitialize
bg_color white
load before.pdb, af3_before
load after.pdb, af3_after
align af3_after, af3_before
hide everything, af3_before or af3_after
show cartoon, af3_before or af3_after
show sticks, (af3_before or af3_after) and not elem H
color red, af3_before
color forest, af3_after
set cartoon_transparency, 0.45, af3_before

# Example: select corrected_before, (af3_before and chain A and resi 4)
# Example: select corrected_after, (af3_after and chain A and resi 4)
select corrected_before, none
select corrected_after, none
show spheres, (corrected_before or corrected_after) and name CA+CB
color orange, corrected_before and name CA+CB
color lime, corrected_after and name CA+CB
label corrected_after and name CA, "corrected %s%s" % (resn, resi)

set sphere_scale, 0.35
set ray_opaque_background, off
group af3_correction, af3_before af3_after corrected_before corrected_after
zoom af3_before or af3_after
