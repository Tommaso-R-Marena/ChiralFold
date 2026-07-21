# Bioinformatics submission package — ChiralFold

**Article type:** Original Paper  
**Journal:** [Bioinformatics](https://academic.oup.com/bioinformatics)  
**Submit at:** https://mc.manuscriptcentral.com/bioinformatics (ScholarOne)

## Contents

| File | Purpose |
|------|---------|
| `chiralfold_bioinformatics.tex` | Main manuscript (OUP template, ~7-page target) |
| `supplementary_material.tex` | Supplementary PDF |
| `references.bib` | Bibliography |
| `oup-authoring-template.cls` | OUP document class |
| `oup-abbrvnat.bst` | Reference style |
| `figures/` | Pre-rendered PNG figures + `generate_figures.py` |
| `data/` | Small CSV/JSON inputs used to regenerate figures |
| `package_overleaf.sh` | Build self-contained Overleaf upload zip |
| `OVERLEAF.md` | Overleaf upload instructions |
| `COVER_LETTER.txt` | Cover letter (customize date & reviewers) |
| `HIGHLIGHTS.md` | Key points for submission form |
| `build_submission.sh` | Compile PDFs and create upload zip |

**Reviewer reproduction path:** [`docs/REVIEWER_COLAB.md`](../../../docs/REVIEWER_COLAB.md) — ordered Colab notebooks + HF Space.

## Build locally

```bash
cd paper/submission/bioinformatics
./build_submission.sh
./check_latex_overflow.sh   # optional: fail if body text spills off page
```

## Overleaf (upload whole directory)

```bash
./package_overleaf.sh   # creates ChiralFold_Bioinformatics_Overleaf.zip
# from the repository root, scripts/build_overleaf_package.sh also copies it to Overleaf/
```

Upload the zip to Overleaf, set main document to `chiralfold_bioinformatics.tex`, and compile. A one-click copy is available at `Overleaf/ChiralFold_Bioinformatics_Overleaf.zip` after running the root helper. All figures are included as PNGs in `figures/`. See `OVERLEAF.md`.

Outputs:
- `chiralfold_bioinformatics.pdf`
- `supplementary_material.pdf`
- `ChiralFold_Bioinformatics_Submission.zip`

## ScholarOne upload checklist

### Manuscript files
- [ ] Main PDF: `chiralfold_bioinformatics.pdf`
- [ ] LaTeX source zip (tex + bib + cls + bst + figures)
- [ ] Supplementary PDF: `supplementary_material.pdf`
- [ ] Cover letter: `COVER_LETTER.txt` → PDF

### Submission form fields
- [ ] **Article type:** Original Paper
- [ ] **Title:** ChiralFold: systematic detection of D-amino acid stereochemistry errors in the Protein Data Bank
- [ ] **Abstract:** copy from `chiralfold_bioinformatics.tex` (structured, ≤150 words)
- [ ] **Keywords:** D-amino acid; chirality; Protein Data Bank; structure validation; AlphaFold 3; stereochemistry
- [ ] **Data availability:** GitHub https://github.com/Tommaso-R-Marena/ChiralFold (v3.5.1); PyPI `chiralfold` (Zenodo DOI optional — see below)

### Zenodo DOI (optional)
A Zenodo DOI is a permanent archive link for a specific GitHub release. Bioinformatics accepts GitHub+PyPI without it.
Skip for initial submission unless your institution requires a DOI. After acceptance (or if you want one now): create a free Zenodo account, connect GitHub, and archive release `v3.5.1`.

### Author information (you must complete in portal)
- [ ] ORCID for Tommaso R. Marena
- [ ] Affiliation: The Catholic University of America, Washington, DC
- [ ] Corresponding author email: marena@cua.edu
- [ ] Funding: None
- [ ] Conflicts: None

### Software availability (Bioinformatics requirement)
- [ ] URL works: https://github.com/Tommaso-R-Marena/ChiralFold
- [ ] PyPI: https://pypi.org/project/chiralfold/
- [ ] Version cited: 3.5.1
- [ ] License: MIT
- [ ] Test data: `results/d_residue_verification.csv`, `tests/fixtures/`
- [ ] Documentation: README.md

### Before you click Submit
- [ ] Confirm main manuscript ≤7 typeset pages (compile and check page count)
- [ ] Verify all figure URLs and GitHub links resolve
- [ ] Replace placeholder suggested reviewers with confirmed names (or remove)
- [ ] Have co-authors approve (if any added)
- [ ] Select open-access option / APC plan if required by your institution

## Page limit note

Bioinformatics Original Papers are limited to **7 printed pages** (~5,000 words excl. figures). The main manuscript is written to be substantive within that cap; extended Lean proofs with derivations, tables, statistics, and figures are in `supplementary_material.tex` (Notes S1). Compile the OUP template and confirm ≤7 pages before submit—do **not** pad past the limit (desk-reject risk above ~20% over). The full unabridged manuscript remains at `paper/chiralfold_paper.tex` for reference.

## After acceptance

- Optional: deposit a Zenodo DOI for GitHub release v3.5.1 (not required if GitHub+PyPI already cited)
- Update manuscript with Bioinformatics DOI
- Add citation to README
