# ChiralFold PyMOL scripts

ChiralFold writes PyMOL `.pml` scripts without importing PyMOL or adding it as
a hard Python dependency.

## Python API

```python
from chiralfold.viz import (
    write_af3_correction_session,
    write_chirality_audit_session,
    write_mirror_comparison_session,
)

write_chirality_audit_session("input.pdb", report, "audit.pml")
write_mirror_comparison_session("l_structure.pdb", "d_structure.pdb", "mirror.pml")
write_af3_correction_session("before.pdb", "after.pdb", "af3_correction.pml")
```

## CLI

```bash
python -m chiralfold.viz.pymol_scripts audit input.pdb report.json -o audit.pml
python -m chiralfold.viz.pymol_scripts mirror l.pdb d.pdb -o mirror.pml
python -m chiralfold.viz.pymol_scripts af3-correction before.pdb after.pdb -o corrected.pml
```

Open interactively with `pymol audit.pml`, or render headlessly with
`pymol -cq audit.pml` after adding a `png output.png, dpi=300` command.

PyMOL is best installed as a system or conda package. The PyPI package is not
portable enough to list as a hard dependency for ChiralFold.
