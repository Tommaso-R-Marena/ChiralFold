# ChiralFold — bioRxiv / arXiv Preprint Submission Package

This directory contains the standalone LaTeX source for the ChiralFold
preprint, formatted for bioRxiv (primary target) and arXiv (secondary).

## Manuscript

- **Title:** ChiralFold: Systematic Detection of D-Amino Acid
  Stereochemistry Errors in the Protein Data Bank
- **Author:** Tommaso R. Marena
- **Institution:** The Catholic University of America, Washington, DC
- **Contact:** marena@cua.edu
- **Main file:** `chiralfold_preprint.tex`
- **Compiled PDF:** `chiralfold_preprint.pdf` (built with TeX Live 2025;
  `pdflatex chiralfold_preprint.tex` × 2 resolves all citations and
  cross-references; no external `.bib` file — the bibliography is inline.)
- **Companion note (separate submission):**
  `../paper/pdb_survey_note.tex` (short PDB survey, ~2,500 words, focused
  on the 12,573-residue audit). Submit as a separate, shorter note if the
  full toolkit manuscript is split.

## Suggested Categories

- **bioRxiv categories:** Biochemistry; Bioinformatics; Synthetic
  Biology. (q-bio.BM as the closest mirror in arXiv terminology.)
- **arXiv categories (if cross-posted):** q-bio.BM (primary); cs.LG
  (secondary, for the AF3 correction component); physics.bio-ph
  (tertiary).
- **Keywords:** D-amino acid; chirality; protein stereochemistry; PDB
  validation; AlphaFold 3; mirror-image phage display; diastereomer;
  drug design.

## Step-by-Step bioRxiv Submission

1. **Recompile from a clean directory.**
   ```bash
   cd preprint
   pdflatex chiralfold_preprint.tex
   pdflatex chiralfold_preprint.tex
   ```
   Confirm `chiralfold_preprint.pdf` has 8 pages, no `??` references, and
   that all figures embed correctly. Required external dirs: only
   `pdb_survey_figures/` (already symlinked / copied in).

2. **Create a bioRxiv account** at https://www.biorxiv.org/submit and
   verify your institutional email (`marena@cua.edu`).

3. **Start a new submission.**
   - Choose **Manuscript Type:** "New Results".
   - Choose **Subject area:** Biochemistry (primary), Bioinformatics
     (secondary).
   - Upload `chiralfold_preprint.pdf` as the main manuscript.
   - Upload the LaTeX source as a single ZIP that contains
     `chiralfold_preprint.tex` and `pdb_survey_figures/`:
     ```bash
     cd preprint
     zip -r chiralfold_preprint_source.zip \
         chiralfold_preprint.tex pdb_survey_figures/
     ```
     bioRxiv stores the source separately from the displayed PDF.

4. **Author and affiliation.** Single author: Tommaso R. Marena,
   The Catholic University of America, Washington, DC, USA. ORCID iD
   (if available) should be added at this step.

5. **Abstract.** Copy the abstract from the LaTeX source verbatim
   (single paragraph, ~260 words). bioRxiv strips LaTeX commands; rerun
   manual cleanup of `$\rho$`, `$\leftrightarrow$`, `$\citep$`, etc., to
   their plain-text equivalents (`rho`, `<->`, "(Childs et al. 2025)").

6. **License.** Recommended: CC BY 4.0. ChiralFold is MIT licensed and
   open source; CC BY is the matching open license for the manuscript.

7. **Competing interests.** Select "No competing interests". Matches
   the statement at the end of the manuscript.

8. **Data Availability statement.** The manuscript's Data Availability
   section already includes GitHub URLs to:
   - `results/d_residue_verification.csv` (12,574 rows)
   - `benchmarks/independent_d_residue_verification.py`
   - `results/error_classification.json`
   - `results/ccd_code_coverage_summary.csv`
   - `results/molprobity_comparison.json`
   These satisfy bioRxiv's data-availability requirement.

9. **Preview, finalise, and submit.** bioRxiv runs an automated
   formatting check; address any flags (typically font embedding or
   page-margin issues). Once approved, you'll receive a DOI of the form
   `10.1101/2026.MM.DD.XXXXXX` within ~48 hours.

10. **Cross-post to arXiv (optional, recommended).** If targeting
    cs.LG / q-bio.BM, the same source zip can be uploaded at
    https://arxiv.org/submit/. arXiv accepts `.tex` directly; no PDF
    upload required. Use category `q-bio.BM` (primary) with cross-list
    to `cs.LG` and `physics.bio-ph`.

## Companion: PDB Survey Note

If the full ChiralFold paper is held back for journal submission, the
shorter PDB-only audit (`../paper/pdb_survey_note.tex`, ~2,500 words,
7 pages) can be submitted to bioRxiv as a standalone note. Suggested
title: "A 12,573-Residue Audit of D-Amino Acid Stereochemistry in the
Protein Data Bank." Same author, same institution, same Data
Availability section.

## Differences from Journal Manuscript

The preprint differs from the journal-submission target only in
metadata: bioRxiv does not require journal-specific class files or
title-page formatting. The LaTeX source here is the standard `article`
class with single-column 1-inch margins — compatible with both
bioRxiv's PDF upload path and arXiv's source compilation.

## Submission Checklist

- [ ] `chiralfold_preprint.pdf` compiles clean (no warnings, no `??`).
- [ ] All figures present in `pdb_survey_figures/`.
- [ ] Abstract under 300 words after LaTeX-to-text conversion.
- [ ] Data Availability links resolve at GitHub (master branch).
- [ ] No author identifiers in figure file metadata.
- [ ] License selected (CC BY 4.0 recommended).
- [ ] Competing interests declared.
- [ ] Single contact email (`marena@cua.edu`).
