# Talk Script: "Finding Hidden Errors in the World's Protein Database"
## Tommaso R. Marena — University Research Day 2026
### 10-Minute Oral Presentation

---

**[SLIDE 1 — Title]** *(30 seconds)*

Good morning, everyone. I'm Tommaso Marena, an undergraduate researcher in the Department of Physics, and today I want to tell you about a detective story — one that starts with a simple question about left and right, and ends with the discovery that some of the most important scientific data in the world has hidden mistakes that nobody noticed for over twenty years.

---

**[SLIDE 2 — Mirror-Image Molecules]** *(60 seconds)*

So let's start with something everyone in this room has — hands. Hold up your hands. Your left hand and your right hand look almost identical, but they're not the same. You can't superimpose them. If you try to put a left glove on your right hand, it doesn't fit.

Molecules have this same property. It's called *chirality*, from the Greek word for "hand." The amino acids that make up every protein in your body come in two mirror-image forms: left-handed and right-handed. And here's the remarkable thing — almost every living thing on Earth uses exclusively left-handed amino acids. That's not a suggestion. It's a rule. Biology is left-handed.

But in the last two decades, pharmaceutical companies have gotten very excited about right-handed amino acids — D-amino acids — because drugs made from them don't get chewed up by your body's enzymes. They last longer, work better, and could treat cancers and infections that resist conventional drugs.

---

**[SLIDE 3 — The Problem]** *(60 seconds)*

Now, when scientists determine the 3D shape of a protein, they deposit it in a global database called the Protein Data Bank — the PDB. It has over 200,000 structures and it's the foundation of modern biology. Every drug company, every research lab, every AI model trained on protein structures uses this database.

There are quality-checking tools that validate these structures before they go in. The most important one is called MolProbity, and it does a great job checking whether left-handed amino acids are built correctly.

But here's the gap: MolProbity was not designed to check right-handed amino acids. D-amino acids. And the best AI protein predictor in the world — AlphaFold 3, from Google DeepMind — gets D-amino acid handedness wrong 51% of the time. That's literally no better than flipping a coin.

So I asked: if nobody is checking, are there errors hiding in the database right now?

---

**[SLIDE 4 — Research Question]** *(30 seconds)*

That's what this project set out to answer. I had three objectives: first, build a tool that can actually check D-amino acid handedness. Second, scan the PDB for errors. And third, figure out whether any errors I found were random or systematic.

---

**[SLIDE 5 — Methods]** *(45 seconds)*

The method turns out to be beautifully simple. Every amino acid has a central carbon atom bonded to four different chemical groups. Those four groups form a tetrahedron — like a little 3D pyramid. And here's the key: the *signed volume* of that pyramid tells you the handedness. If the volume is negative, it's right-handed — D. If it's positive, it's left-handed — L.

So for any residue labeled as a D-amino acid in the PDB, I just compute this signed volume from the raw atomic coordinates. If I get a positive number, that means the coordinates say "left-handed" even though the label says "right-handed." That's an error.

No AI. No machine learning. Just basic geometry applied to every D-amino acid in the database.

---

**[SLIDE 6 — 21 Errors Found]** *(45 seconds)*

Here's what I found. Out of 1,677 D-amino acid residues across 231 PDB files, 21 have the wrong handedness. They're labeled D, but their coordinates are actually L. That's a 1.3% error rate across 12 distinct structures.

Now, 1.3% might sound small. But remember — these are in the foundational database of structural biology. Every computational study that uses these structures is using wrong stereochemistry for those residues. And these errors are invisible to existing tools.

---

**[SLIDE 7 — All Before 2006]** *(60 seconds)*

But here's where it gets really interesting. When I looked at *when* these structures were deposited, a pattern emerged. Every single one of the 12 error-containing structures was deposited between 1992 and 2005. Zero errors in anything deposited after 2005. The statistical test gives a p-value of 0.003 — this is not random.

What happened in 2006? The worldwide PDB organization launched a massive "remediation" effort — essentially a cleanup of the entire archive. They fixed many stereochemistry issues. But they didn't specifically check D-amino acid labels. So these errors slipped through, and they've been sitting there for twenty years.

---

**[SLIDE 8 — The 1HHZ Case]** *(45 seconds)*

Let me show you the most striking example. PDB structure 1HHZ was determined at 0.99 Angstrom resolution. That is atomic resolution — you can literally see individual atoms in the electron density map. The data quality is essentially perfect.

But the D-Alanine at position 1 in chain E? Its coordinates clearly show left-handed stereochemistry. Signed volume: positive 2.70. The data is flawless. The label is wrong. This single case proves that these errors are not about bad data or low resolution — they're about mislabeling during deposition.

---

**[SLIDE 9 — Mirror-Image Drug Design]** *(60 seconds)*

Beyond finding errors, ChiralFold also enables something practical: mirror-image drug design. The idea is simple. If you have a left-handed peptide that binds to a cancer protein, you can mirror it to create a right-handed version that binds the same target — but now it's resistant to the body's enzymes.

I demonstrated this on MDM2, a protein that suppresses the tumor suppressor p53. I took the crystal structure of p53 bound to MDM2, mirrored the peptide, and showed that the three critical binding residues — Phenylalanine, Tryptophan, and Leucine — are perfectly preserved as D-amino acids. The mirror transformation maintains 105 atomic contacts and 10 hydrogen bonds with zero error. And this matches the geometry of a real D-peptide drug — dPMI-gamma — that binds MDM2 at 53 nanomolar affinity.

---

**[SLIDE 10 — AlphaFold 3]** *(45 seconds)*

I also benchmarked against AlphaFold 3 using the actual D-peptide sequences from a 2025 study by Childs and colleagues. These are real sequences from real crystal structures. AlphaFold 3 produces chirality errors 51% of the time on these sequences. ChiralFold: zero percent. And ChiralFold can automatically detect and correct the errors in AlphaFold's outputs.

This isn't because ChiralFold is "smarter" — it's because it builds chirality in from the start, by construction, rather than hoping a neural network learns it from data.

---

**[SLIDE 11 — What This Means]** *(45 seconds)*

So what does this all mean? Three things.

First, there are genuine, verifiable stereochemistry errors hiding in the Protein Data Bank from the pre-remediation era. These are not theoretical — I can give you twelve PDB IDs right now that you can check yourself.

Second, existing validation tools have a blind spot for D-amino acids. ChiralFold fills that gap as an open-source, pip-installable Python package that anyone can use.

And third, as D-peptide drugs become more important — and this is now a multi-billion dollar field — having correct stereochemistry in our structural databases isn't optional. It's essential.

---

**[SLIDE 12 — Conclusions]** *(30 seconds)*

To summarize: I built ChiralFold, the first tool specifically designed to validate D-amino acid stereochemistry. Using it, I found 21 genuine errors in 12 PDB structures — all from before the 2006 remediation, invisible to MolProbity, spanning 0.99 to 2.70 Angstrom resolution. The tool also enables mirror-image drug design and can correct AlphaFold 3 outputs. It's freely available on GitHub.

Thank you. I'm happy to take questions.

---

**Total estimated time: 9 minutes 30 seconds**
*(At ~150 words/minute speaking pace, this script is approximately 1,400 words)*
