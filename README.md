# Structural Biology: Coursework, Tools & Semester Project

Work from a semester course in Structural Biology (BS–MS, IISER Bhopal, instructed by Professor Sunando Datta): a
backbone-geometry visualiser built from scratch, PyMOL code guides written while
learning the program, and a
computational semester project on whether AlphaFold2 can see the active-site
change that reprograms an enzyme's product.



---

## Contents

```
.
├── torsion-visualizer/       Interactive phi/psi/omega visualiser (single-file HTML)
├── pymol/                    guided scripts for learning basics of PyMol          
└── eizs-alphafold-project/   Semester project: EIZS active-site mutants vs AlphaFold2
    ├── eizs_analysis.py
    ├── report/               compiled PDF
    ├── protocol/             Step-by-step ColabFold execution protocol
    └── results/
```



---

## 1. Protein Backbone Torsion Angle Visualiser

`torsion-visualizer/`

A single-file HTML tool for building a poly-alanine backbone residue by residue
and watching the geometry respond to φ, ψ and ω. You can directly open the file in a browser.

**Relevance.** Ramachandran plots and Newman projections are usually taught
as static pictures. Here they are driven by actual 3-D coordinates, so the
relationship between a dihedral slider, the Newman projection down that bond,
and the position on the Ramachandran map is visible at once.

**What it does**

- 3-D canvas viewport with drag-rotate, zoom, pan and depth shading, implemented
  without WebGL
- Atom placement via **NeRF** (Natural Extension Reference Frame) using
  **Engh & Huber** ideal bond lengths and angles
- Per-residue φ/ψ/ω sliders; chain mode up to 10 residues with "link all" to
  propagate one conformation along the chain
- Live Newman projections down **N→Cα** and **Cα→C′**, projected from the real
  coordinates rather than drawn schematically
- Ramachandran map with favoured/allowed overlays and **drag-to-set** — click a
  region and the 3-D structure adopts it
- Toggles: atoms, labels, hydrogens, Cβ, terminal caps (acetyl / N-methyl),
  amide planes, φ/ψ bond highlighting with torsion arcs, steric clash detection
  on van der Waals radii
- **PDB export** of the current conformation

**Secondary-structure regions covered:** α-helix, β-strand, PPII, 3₁₀,
left-handed α.

**Note:** I do not know how to code this in HTML :), this was made using Claude.

---

## 2. Guided Scripts for PyMol

`pymol/`

Written while learning PyMOL from scratch for the practical component of the
course. 

**Tasks covered:** structure loading and cleaning · selections and colouring ·
distances, including salt bridges and π-stacking · bond angles · torsion angles
(φ/ψ/ω/χ₁) · protein–ligand interactions with hydrogen-bond listing ·
protein–peptide interface analysis · superposition and RMSD · surface area and
SASA · publication figure generation · buried surface area (BSA) for both
protein–peptide and protein–small-molecule cases, by the FreeSASA route and a
native-PyMOL route.

**Structures used:** `1UBQ` (ubiquitin) · `1STP` (streptavidin–biotin) · `1YCR`
(MDM2–p53 peptide) · `2B2W`.


---

## 4. Semester Project — Can AlphaFold2 See a Reprogrammed Active Site?

`eizs-alphafold-project/`

### The question

***epi*-isozizaene synthase (EIZS)** from *Streptomyces coelicolor* A3(2) takes
one achiral substrate (farnesyl diphosphate) and makes one tricyclic product.
Single substitutions in its active site turn it into six *different*
sesquiterpene synthases:

| Variant | Becomes |
|---|---|
| F96A | (*E*)-β-farnesene synthase |
| F96W | zizaene synthase |
| F95H | β-curcumene synthase |
| F95M | β-acoradiene synthase |
| F198L | β-cedrene synthase |
| F96V, W203F | (*Z*)-γ-bisabolene synthase |

None of these residues is catalytic. They line a hydrophobic contour (Y69, F95,
F96, F198, W203) that acts as a template for a carbocation cascade — change the
shape of the template, change the product.

AlphaFold2 sees only sequence. So: **does a single point substitution register
in the predicted structure at all, and if it does, does it register in a way
that tracks the change in product?**


## Notes

Study guides and PyMOL material are personal course notes, not official course
documents. The project report is coursework and is included as submitted. Code
in this repository is free to reuse; please cite the papers above rather than
this repository for any of the underlying science.