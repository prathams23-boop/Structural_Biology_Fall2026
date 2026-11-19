#!/usr/bin/env python3
"""Rebuild a pLDDT plot from saved ColabFold models.

ColabFold writes each residue's pLDDT into the B-factor column of the output
PDB, so the confidence trace can be recovered from the models alone, with no
need to re-run anything.

Usage:
    python plot_plddt.py models/            # folder of <VARIANT>.pdb
    python plot_plddt.py models/ --out plddt.png
"""
import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ACTIVE_SITE = [69, 95, 96, 198, 203, 325]


def plddt_from_pdb(path):
    """Per-residue pLDDT, read off the CA atom of each residue."""
    resnum, plddt = [], []
    seen = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue
            try:
                n = int(line[22:26])
                b = float(line[60:66])
            except ValueError:
                continue
            if n in seen:
                continue
            seen.add(n)
            resnum.append(n)
            plddt.append(b)
    return np.array(resnum), np.array(plddt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("models", help="folder containing <VARIANT>.pdb")
    ap.add_argument("--out", default="plddt_all_variants.png")
    ap.add_argument("--site", default=",".join(map(str, ACTIVE_SITE)))
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.models, "*.pdb")))
    if not paths:
        raise SystemExit(f"no PDB files in {args.models}")
    site = [int(x) for x in args.site.split(",") if x.strip()]

    fig, ax = plt.subplots(figsize=(9, 4))
    for p in paths:
        name = os.path.splitext(os.path.basename(p))[0]
        n, b = plddt_from_pdb(p)
        if len(n) == 0:
            print(f"  {name}: no CA atoms found, skipped")
            continue
        ax.plot(n, b, lw=1.3, label=f"{name}  (mean {b.mean():.1f})")
        print(f"  {name:<8} mean pLDDT {b.mean():6.2f}   min {b.min():6.2f}   "
              f"active site {np.mean([b[n == r][0] for r in site if (n == r).any()]):6.2f}")

    for r in site:
        ax.axvline(r, color="0.75", ls=":", lw=0.8, zorder=0)
    ax.axhline(90, color="g", ls="--", lw=0.7)
    ax.axhline(70, color="orange", ls="--", lw=0.7)
    ax.set_xlabel("Residue number")
    ax.set_ylabel("pLDDT")
    ax.set_ylim(0, 100)
    ax.set_title("AlphaFold2 confidence per residue; dotted lines mark active-site residues")
    ax.legend(fontsize=8, loc="lower center", ncol=3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()