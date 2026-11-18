# Does AlphaFold2 reproduce the real active-site change?

Reference frame and ligand source: **WT.pdb** MG (3 atoms), POP (9 atoms), BTM (14 atoms). Every structure is superposed onto it and shares this one ligand set, so all volume differences are differences in protein atoms.


## Table 1 - Every structure

| Structure | Residues | CA RMSD to ref (A) | Pocket volume (A^3) | Min ligand-protein distance (A) | Contact residues |
|---|---|---|---|---|---|
| AF2:WT | 361 | 0.72 | 214 | 1.99 **CLASH** | 22 |
| XTAL:WT | 340 | 0.00 | 141 | 2.06 **CLASH** | 24 |
| AF2:A236G | 361 | 0.70 | 211 | 1.91 **CLASH** | 22 |
| XTAL:A236G | 340 | 0.16 | 143 | 1.98 **CLASH** | 24 |
| AF2:F95H | 361 | 0.61 | 187 | 1.94 **CLASH** | 22 |
| XTAL:F95H | 340 | 0.22 | 134 | 2.07 **CLASH** | 24 |
| AF2:F95M | 361 | 0.61 | 226 | 1.97 **CLASH** | 22 |
| XTAL:F95M | 340 | 0.17 | 149 | 2.01 **CLASH** | 24 |
| AF2:L72V | 361 | 0.68 | 214 | 1.97 **CLASH** | 22 |
| XTAL:L72V | 340 | 0.13 | 147 | 1.90 **CLASH** | 23 |

Any value below 2.2 A in the fifth column means the transplanted ligand overlaps protein in that structure and its pocket volume must not be used.


## Table 2 - THE RESULT: predicted vs real change

| Mutant | dV side chain (A^3) | d Pocket, experiment (A^3) | d Pocket, AlphaFold2 (A^3) | Residual (A^3) | Same sign? |
|---|---|---|---|---|---|
| A236G | -28.5 | +2.4 | -3.9 | -6.3 | **no** |
| F95H | -36.7 | -7.3 | -27.9 | -20.6 | yes |
| F95M | -27.0 | +8.0 | +11.0 | +3.0 | yes |
| L72V | -26.7 | +6.3 | -0.2 | -6.5 | **no** |

*d Pocket, experiment* = mutant crystal minus WT crystal: what the mutation really does. *d Pocket, AlphaFold2* = mutant model minus WT model: what AlphaFold2 predicts from sequence alone. *Residual* is the gap between them. *dV side chain* is the null model -- the volume the substitution adds or removes with nothing else moving.


## Table 3 - Minimum residue-ligand distance (A)

| Structure | 69 | 95 | 96 | 198 | 203 | 325 |
|---|---|---|---|---|---|---|
| AF2:WT | Y69: 7.10 | F95: 2.62 | F96: 3.48 | F198: 2.83 | W203: 3.76 | W325: 3.78 |
| XTAL:WT | Y69: 6.69 | F95: 3.57 | F96: 3.65 | F198: 3.16 | W203: 4.06 | W325: 3.84 |
| AF2:A236G | Y69: 7.09 | F95: 2.64 | F96: 3.51 | F198: 2.87 | W203: 3.82 | W325: 3.61 |
| XTAL:A236G | Y69: 6.72 | F95: 3.51 | F96: 3.69 | F198: 3.11 | W203: 4.02 | W325: 3.85 |
| AF2:F95H | Y69: 7.01 | H95: 2.52 | F96: 3.37 | F198: 2.80 | W203: 3.80 | W325: 4.00 |
| XTAL:F95H | Y69: 6.77 | H95: 3.47 | F96: 3.79 | F198: 3.08 | W203: 3.96 | W325: 3.79 |
| AF2:F95M | Y69: 7.02 | M95: 3.47 | F96: 3.63 | F198: 2.81 | W203: 3.86 | W325: 4.31 |
| XTAL:F95M | Y69: 6.66 | M95: 3.36 | F96: 3.73 | F198: 3.20 | W203: 3.95 | W325: 3.74 |
| AF2:L72V | Y69: 7.18 | F95: 2.64 | F96: 3.60 | F198: 2.84 | W203: 3.75 | W325: 3.95 |
| XTAL:L72V | Y69: 6.67 | F95: 3.43 | F96: 3.76 | F198: 3.05 | W203: 3.99 | W325: 3.89 |

Check the XTAL:WT row shows F95, F96, F198, W203. If it does not, the crystal numbering does not match and nothing above is valid.


## Table 4 - Residues contacting the ligands, WT crystal

| Residue | Contacts | Min distance (A) |
|---|---|---|
| F95 | 10 | 3.57 |
| F198 | 8 | 3.16 |
| F96 | 8 | 3.65 |
| N240 | 6 | 2.07 |
| R338 | 5 | 2.68 |
| E248 | 4 | 2.06 |
| D99 | 4 | 2.09 |
| Y339 | 3 | 2.56 |
| R194 | 3 | 2.75 |
| Y172 | 3 | 3.85 |
| S244 | 2 | 2.32 |
| K247 | 2 | 2.88 |
| T197 | 2 | 3.43 |
| H333 | 2 | 3.56 |
| F332 | 2 | 3.62 |
| L72 | 2 | 3.68 |
| W325 | 2 | 3.84 |
| E175 | 2 | 4.03 |
| W203 | 2 | 4.06 |
| A199 | 2 | 4.26 |
