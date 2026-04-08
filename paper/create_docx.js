const fs = require('fs');
const docx = require('docx');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Header, Footer, AlignmentType, HeadingLevel, BorderStyle, WidthType,
  ShadingType, PageNumber, ExternalHyperlink,
} = docx;

// Borders
const thinBorder = { style: BorderStyle.SINGLE, size: 1, color: '999999' };
const allBorders = { top: thinBorder, bottom: thinBorder, left: thinBorder, right: thinBorder };

// Helper functions
function bodyPara(runs, opts = {}) {
  return new Paragraph({
    spacing: { after: 200, line: 276 },
    ...opts,
    children: runs,
  });
}

function bodyText(text, opts = {}) {
  return new TextRun({ text, font: 'Times New Roman', size: 22, ...opts }); // 11pt = 22 half-pts
}

function heading1(text) {
  return new Paragraph({
    spacing: { before: 360, after: 200 },
    children: [new TextRun({ text, font: 'Times New Roman', size: 28, bold: true })], // 14pt
  });
}

function heading2(text) {
  return new Paragraph({
    spacing: { before: 280, after: 160 },
    children: [new TextRun({ text, font: 'Times New Roman', size: 24, bold: true })], // 12pt
  });
}

// Table 1 data
const tableData = [
  ['PDB ID', 'Residue', 'Chain', 'Position', 'Signed Vol.', 'Resolution', 'Dep. Date', 'Errors'],
  ['1ABI', 'DPN', 'I', '56', '+2.49', '2.30 \u00C5', '1992', '1'],
  ['1BG0', 'DAR', 'A', '403', '+2.58', '1.86 \u00C5', '1998', '1'],
  ['1D7T', 'DTY', 'A', '4', '+1.85', 'NMR', '1999', '1'],
  ['1HHZ', 'DAL', 'E', '1', '+2.70', '0.99 \u00C5', '2000', '1'],
  ['1KO0', 'DLY', 'A', '542', '+0.12', '2.20 \u00C5', '2001', '1'],
  ['1MCB', 'DHI', 'P', '3', '+2.60', '2.70 \u00C5', '1993', '1'],
  ['1OF6', 'DTY', 'A\u2013H', '1369\u20131370', '+2.51 to +2.67', '2.10 \u00C5', '2003', '8'],
  ['1P52', 'DAR', 'A', '403', '+2.54', '1.90 \u00C5', '2003', '1'],
  ['1UHG', 'DSN', 'D', '164', '+2.21', '1.90 \u00C5', '2003', '1'],
  ['1XT7', 'DSG', 'A', '3', '+2.55', 'NMR', '2004', '1'],
  ['2AOU', 'DCY', 'A', '248', '+2.67', '2.30 \u00C5', '2005', '1'],
  ['2ATS', 'DLY', 'A', '3001\u20133003', '+2.56 to +2.59', '1.90 \u00C5', '2005', '3'],
];

const colWidths = [1000, 1000, 800, 1200, 1600, 1200, 1100, 860];
const tableWidth = colWidths.reduce((a, b) => a + b, 0);

function makeCell(text, isHeader, colIdx) {
  const cellChildren = [
    new Paragraph({
      spacing: { after: 0 },
      children: [
        new TextRun({
          text,
          font: 'Times New Roman',
          size: 18, // 9pt for table
          bold: isHeader,
        }),
      ],
    }),
  ];
  return new TableCell({
    borders: allBorders,
    width: { size: colWidths[colIdx], type: WidthType.DXA },
    shading: isHeader ? { fill: 'D9E2F3', type: ShadingType.CLEAR } : undefined,
    margins: { top: 40, bottom: 40, left: 60, right: 60 },
    children: cellChildren,
  });
}

const tableRows = tableData.map((row, rowIdx) =>
  new TableRow({
    children: row.map((cell, colIdx) => makeCell(cell, rowIdx === 0, colIdx)),
  })
);

const table1 = new Table({
  width: { size: tableWidth, type: WidthType.DXA },
  columnWidths: colWidths,
  rows: tableRows,
});

// Build document
const doc = new Document({
  styles: {
    default: {
      document: {
        run: { font: 'Times New Roman', size: 22 },
      },
    },
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
        },
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun({ children: [PageNumber.CURRENT], font: 'Times New Roman', size: 20 }),
              ],
            }),
          ],
        }),
      },
      children: [
        // Title
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 120 },
          children: [
            new TextRun({
              text: 'ChiralFold: Systematic Detection of D-Amino Acid Stereochemistry Errors in the Protein Data Bank',
              font: 'Times New Roman',
              size: 36, // 18pt
              bold: true,
            }),
          ],
        }),

        // Author
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 40 },
          children: [
            new TextRun({
              text: 'Tommaso R. Marena',
              font: 'Times New Roman',
              size: 24, // 12pt
              italics: true,
            }),
          ],
        }),

        // Affiliation
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 40 },
          children: [
            new TextRun({
              text: 'The Catholic University of America, Washington, DC',
              font: 'Times New Roman',
              size: 22,
              italics: true,
            }),
          ],
        }),

        // Email
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 300 },
          children: [
            new TextRun({
              text: 'marena@cua.edu',
              font: 'Times New Roman',
              size: 22,
            }),
          ],
        }),

        // Abstract heading
        new Paragraph({
          spacing: { after: 120 },
          children: [
            new TextRun({ text: 'Abstract', font: 'Times New Roman', size: 24, bold: true }),
          ],
        }),

        // Abstract text
        bodyPara([
          bodyText('We present ChiralFold, a general-purpose protein stereochemistry toolkit that provides chirality-correct coordinate generation, PDB structure auditing, and mirror-image L\u2194D transformation. An independent verification of 1,677 D-amino acid residues across 231 PDB files \u2014 using only the signed tetrahedron volume at each C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' \u2014 identified 21 genuine chirality errors (1.3% error rate) in 12 deposited structures where C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' coordinates show L-stereochemistry despite D-amino acid labels. All 12 error-containing structures were deposited between 1992 and 2005, with zero errors in post-2005 entries (Mann\u2013Whitney U=278, p=0.0027), consistent with the 2006\u20132008 wwPDB remediation. Resolution does not predict errors (p=0.19); errors span 0.99\u20132.70 \u00C5, indicating a labeling problem rather than a data quality problem. ChiralFold\u2019s auditor, calibrated against wwPDB/MolProbity validation reports on 31 structures (Ramachandran outlier Spearman \u03C1=0.49, p=0.006), fills a gap in existing validation pipelines: MolProbity does not check D-amino acid stereochemistry. The toolkit also includes an AlphaFold 3 chirality correction pipeline addressing the documented 51% violation rate on D-peptides (Childs et al., 2025), and a mirror-image binder design workflow validated on the MDM2:dPMI-\u03B3 system (K'),
          bodyText('d', { subScript: true }),
          bodyText('=53 nM). ChiralFold is available as a pip-installable Python package at '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'https://github.com/Tommaso-R-Marena/ChiralFold', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://github.com/Tommaso-R-Marena/ChiralFold',
          }),
          bodyText('.'),
        ]),

        // 1. Introduction
        heading1('1. Introduction'),

        bodyPara([
          bodyText('D-amino acids play an increasingly important role in drug design due to their resistance to proteolysis and improved bioavailability. Mirror-image phage display has produced D-peptide therapeutics with nanomolar affinity, including dPMI-\u03B3 (K'),
          bodyText('d', { subScript: true }),
          bodyText('=53 nM against MDM2). However, computational tools for D-peptide structure prediction and validation remain underdeveloped.'),
        ]),

        bodyPara([
          bodyText('AlphaFold 3, the current state of the art for protein structure prediction, produces a 51% per-residue chirality violation rate on D-peptides \u2014 equivalent to random assignment (Childs, Zhou & Donald, 2025). This fundamental limitation arises from AF3\u2019s diffusion architecture, which denoises atom coordinates without enforcing hard stereochemical constraints.'),
        ]),

        bodyPara([
          bodyText('Existing validation tools also have a blind spot. MolProbity, the standard for PDB structure quality assessment, validates C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' chirality for standard L-amino acids but does not specifically check D-amino acid stereochemistry. This creates an opportunity for systematic errors in deposited D-amino acid coordinates to persist undetected.'),
        ]),

        bodyPara([
          bodyText('We present ChiralFold, a pip-installable Python toolkit that addresses both gaps: it provides guaranteed chirality-correct D-peptide coordinate generation and a comprehensive PDB auditor that validates D-amino acid stereochemistry. Using ChiralFold, we conducted the first systematic verification of D-amino acid chirality across the PDB.'),
        ]),

        // 2. Methods
        heading1('2. Methods'),

        heading2('2.1 Signed Tetrahedron Volume'),

        bodyPara([
          bodyText('Chirality at each C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' is determined by the signed volume of the tetrahedron formed by the four substituents:'),
        ]),

        // Equation
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { before: 200, after: 200 },
          children: [
            bodyText('signed_volume = (N \u2212 C'),
            bodyText('\u03B1', { subScript: true }),
            bodyText(') \u00B7 ((C \u2212 C'),
            bodyText('\u03B1', { subScript: true }),
            bodyText(') \u00D7 (C'),
            bodyText('\u03B2', { subScript: true }),
            bodyText(' \u2212 C'),
            bodyText('\u03B1', { subScript: true }),
            bodyText('))'),
          ],
        }),

        bodyPara([
          bodyText('For L-amino acids (S-configuration), this volume is negative by CIP convention. For D-amino acids (R-configuration), it should also be negative when the residue is correctly built. A positive signed volume at a D-labeled residue indicates that the deposited coordinates have L-stereochemistry \u2014 a chirality error.'),
        ]),

        bodyPara([
          bodyText('This method requires only the four backbone atom positions (N, C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(', C, C'),
          bodyText('\u03B2', { subScript: true }),
          bodyText(') and uses no force field or external dependencies beyond numpy.'),
        ]),

        heading2('2.2 Independent Verification Protocol'),

        bodyPara([
          bodyText('To ensure our findings are not artifacts of ChiralFold\u2019s implementation, the D-residue verification was performed independently: the script uses only numpy for the cross product and dot product, parses PDB files with standard string operations, and imports no ChiralFold code. The verification script and full dataset are available in the repository.'),
        ]),

        heading2('2.3 PDB Auditor'),

        bodyPara([
          bodyText('ChiralFold\u2019s auditor validates six quality criteria: C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' chirality (signed volume), bond geometry (RMSD to Engh & Huber ideal values), Ramachandran backbone dihedrals (hybrid empirical PDB grid + calibrated rectangular regions), peptide bond planarity (\u03C9 deviation from 180\u00B0), steric clashes (van der Waals overlap with backbone hydrogen placement), and a weighted composite score (0\u2013100).'),
        ]),

        heading2('2.4 Mirror-Image Pipeline'),

        bodyPara([
          bodyText('The L\u2194D transformation reflects all atomic coordinates across a plane: (x,y,z) \u2192 (\u2212x,y,z). This inverts all stereocenters (det(R) = \u22121), preserving bond lengths, bond angles, and torsion magnitudes exactly. Residue names are mapped to their D-amino acid equivalents (ALA\u2192DAL, TRP\u2192DTR, etc.).'),
        ]),

        // 3. Results
        heading1('3. Results'),

        heading2('3.1 D-Amino Acid Chirality Errors in the PDB'),

        bodyPara([
          bodyText('We verified 1,677 D-amino acid residues with complete backbone atoms (N, C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(', C, C'),
          bodyText('\u03B2', { subScript: true }),
          bodyText(') across 231 PDB files. Of these, 1,656 (98.7%) showed correct D-chirality (negative signed volume) and 21 (1.3%) showed L-chirality (positive signed volume) despite being labeled as D-amino acids.'),
        ]),

        bodyPara([
          bodyText('The 21 errors occur in 12 PDB structures:'),
        ]),

        // Table 1
        table1,

        // Table caption
        new Paragraph({
          spacing: { before: 80, after: 240 },
          children: [
            new TextRun({
              text: 'Table 1. The 12 PDB structures containing D-amino acid chirality errors. Signed volume is in \u00C5\u00B3; positive values indicate L-stereochemistry at the C',
              font: 'Times New Roman',
              size: 18,
              italics: true,
            }),
            new TextRun({
              text: '\u03B1',
              font: 'Times New Roman',
              size: 18,
              italics: true,
              subScript: true,
            }),
            new TextRun({
              text: '.',
              font: 'Times New Roman',
              size: 18,
              italics: true,
            }),
          ],
        }),

        heading2('3.2 Errors Correlate with Pre-Remediation Deposition'),

        bodyPara([
          bodyText('All 12 error-containing structures were deposited between 1992 and 2005. Zero errors were found in structures deposited after 2005. Deposition year significantly predicts the presence of errors (Mann\u2013Whitney U=278, p=0.0027). This temporal pattern is consistent with the 2006\u20132008 wwPDB remediation effort, which systematically corrected stereochemistry assignments across the archive but did not specifically address D-amino acid residue labeling.'),
        ]),

        bodyPara([
          bodyText('Resolution does not predict errors (Mann\u2013Whitney U=112, p=0.19). The highest-resolution error structure is 1HHZ at 0.99 \u00C5 \u2014 an ultra-high-resolution structure where the electron density is unambiguous but the residue was labeled as D-Alanine (DAL) while the C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' coordinates show L-stereochemistry (signed volume +2.70). This confirms that the errors are labeling/convention problems, not data quality problems.'),
        ]),

        heading2('3.3 MolProbity Comparison'),

        bodyPara([
          bodyText('ChiralFold\u2019s auditor was benchmarked against wwPDB/MolProbity validation reports on 31 PDB structures spanning 0.48\u20133.4 \u00C5 resolution (X-ray, NMR, cryo-EM). Ramachandran outlier percentages showed significant agreement (Spearman \u03C1=0.49, p=0.006), with ChiralFold reporting 0.60% mean outliers vs wwPDB\u2019s 0.64%.'),
        ]),

        bodyPara([
          bodyText('MolProbity is stronger in Ramachandran contour precision (data-derived from ~100K structures vs ChiralFold\u2019s calibrated rectangles + empirical grid) and rotamer completeness. ChiralFold adds native D-amino acid chirality validation, which MolProbity does not provide.'),
        ]),

        heading2('3.4 AlphaFold 3 Chirality on Real D-Peptide Sequences'),

        bodyPara([
          bodyText('Using the actual D-peptide sequences from Childs et al. (2025) \u2014 DEHELLETAARWFYEIAKR (PDB 7YH8, DP19), LWQHEATWK (PDB 5N8T, DP9), and DWWPLAFEALLR (PDB 3IWY, DP12) \u2014 ChiralFold produces 0/478 chirality violations across 41 sequence/pattern combinations. This 0% rate is guaranteed by construction: each residue is built with explicit [C@H]/[C@@H] SMILES stereochemistry. AF3\u2019s 51% rate reflects its diffusion model treating D-residues as noise.'),
        ]),

        heading2('3.5 Mirror-Image Binder Design'),

        bodyPara([
          bodyText('The p53:MDM2 crystal structure (PDB 1YCR) was mirrored to generate a D-peptide binder candidate. The binding triad (Phe19/Trp23/Leu26) is preserved as D-Phe/D-Trp/D-Leu \u2014 the same hotspot residues that the experimental dPMI-\u03B3 (K'),
          bodyText('d', { subScript: true }),
          bodyText('=53 nM, PDB 3IWY) uses. All 13 backbone \u03C6 angles are exactly sign-inverted. Contact geometry is preserved by construction: 105 C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText('\u2013C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' contacts and 10 H-bond donor\u2013acceptor pairs maintained within 0.001 \u00C5 across the L\u2192D transformation.'),
        ]),

        // 4. Discussion
        heading1('4. Discussion'),

        bodyPara([
          bodyText('The 1.3% error rate in deposited D-amino acid stereochemistry is small but systematic. The temporal clustering before 2006 suggests these errors were introduced during an era when D-amino acid validation was not part of standard deposition pipelines. The wwPDB remediation addressed many stereochemistry issues but these D-residue errors persisted because no existing tool specifically validated D-amino acid C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' chirality.'),
        ]),

        bodyPara([
          bodyText('The 1HHZ case is particularly illustrative: at 0.99 \u00C5 resolution, the electron density is unambiguous, yet the residue is mislabeled. This rules out resolution as a confound and points to a convention/software error during deposition.'),
        ]),

        bodyPara([
          bodyText('ChiralFold fills this validation gap as a lightweight, pip-installable Python tool. While it does not replace MolProbity for comprehensive structure validation, it provides the only available programmatic D-amino acid chirality check and a mirror-image coordinate generation pipeline relevant to the growing field of D-peptide drug design.'),
        ]),

        // 5. Data Availability
        heading1('5. Data Availability'),

        bodyPara([
          bodyText('The complete verification dataset (1,678 rows with raw N/C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText('/C/C'),
          bodyText('\u03B2', { subScript: true }),
          bodyText(' coordinates for every D-amino acid residue) is available at: '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'https://github.com/Tommaso-R-Marena/ChiralFold/blob/master/results/d_residue_verification.csv', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://github.com/Tommaso-R-Marena/ChiralFold/blob/master/results/d_residue_verification.csv',
          }),
        ]),

        bodyPara([
          bodyText('The independent verification script (no ChiralFold dependencies) is at: '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'https://github.com/Tommaso-R-Marena/ChiralFold/blob/master/benchmarks/independent_d_residue_verification.py', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://github.com/Tommaso-R-Marena/ChiralFold/blob/master/benchmarks/independent_d_residue_verification.py',
          }),
        ]),

        // References
        heading1('References'),

        bodyPara([
          bodyText('1. Childs CM, Zhou P, Donald BR. Has AlphaFold 3 Solved the Protein Folding Problem for D-Peptides? '),
          bodyText('bioRxiv', { italics: true }),
          bodyText(' 2025.03.14.643307 (2025). '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'doi:10.1101/2025.03.14.643307', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://doi.org/10.1101/2025.03.14.643307',
          }),
        ]),

        bodyPara([
          bodyText('2. Liu M et al. D-peptide inhibitors of the p53\u2013MDM2 interaction for targeted molecular therapy of malignant neoplasms. '),
          bodyText('Proc Natl Acad Sci', { italics: true }),
          bodyText(' 107:14321\u201314326 (2010). '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'doi:10.1073/pnas.1008930107', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://doi.org/10.1073/pnas.1008930107',
          }),
        ]),

        bodyPara([
          bodyText('3. Chen VB et al. MolProbity: all-atom structure validation for macromolecular crystallography. '),
          bodyText('Acta Cryst', { italics: true }),
          bodyText(' D66:12\u201321 (2010). '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'doi:10.1107/S0907444909042073', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://doi.org/10.1107/S0907444909042073',
          }),
        ]),

        bodyPara([
          bodyText('4. Lovell SC et al. Structure validation by C'),
          bodyText('\u03B1', { subScript: true }),
          bodyText(' geometry: \u03C6,\u03C8 and C'),
          bodyText('\u03B2', { subScript: true }),
          bodyText(' deviation. '),
          bodyText('Proteins', { italics: true }),
          bodyText(' 50:437\u2013450 (2003). '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'doi:10.1002/prot.10286', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://doi.org/10.1002/prot.10286',
          }),
        ]),

        bodyPara([
          bodyText('5. Henrick K et al. Remediation of the protein data bank archive. '),
          bodyText('Nucleic Acids Res', { italics: true }),
          bodyText(' 36:D426\u2013D433 (2008). '),
          new ExternalHyperlink({
            children: [new TextRun({ text: 'doi:10.1093/nar/gkm937', style: 'Hyperlink', font: 'Times New Roman', size: 22 })],
            link: 'https://doi.org/10.1093/nar/gkm937',
          }),
        ]),
      ],
    },
  ],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync('/home/user/workspace/chiralfold/paper/chiralfold_paper.docx', buf);
  console.log('DOCX created successfully');
});
