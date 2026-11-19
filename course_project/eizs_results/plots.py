#!/usr/bin/env python3
"""Extra figures from the models and crystals you already have.

Makes two plots, neither of which needs Colab, Biopython or a re-run:

  dev_per_residue.png  where along the chain each AlphaFold2 model deviates
                       from its own crystal structure, after superposition
  site_error.png       how far off the model is for each active-site residue,
                       measured as distance to the bound ligand

Usage, from inside the eizs_results folder:
    python3 plots.py
or with explicit paths:
    python3 plots.py --models models --crystals crystals --results results
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SITE = [69, 95, 96, 198, 203, 325]
VARIANTS = ["WT", "F95H", "F95M", "L72V", "A236G"]


def read_ca(path):
    """{resseq: (resname, xyz)} for every CA atom."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA" or line[16] not in (" ", "A"):
                continue
            try:
                n = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            out.setdefault(n, (line[17:20].strip(), xyz))
    return out


def best_offset(a, b):
    """Numbering offset that makes the two sets agree on the most residues.

    Their numbering already matches, so this returns 0, but scanning costs
    nothing and keeps the script working if a file is numbered differently.
    """
    best, best_n = 0, -1
    for off in range(-60, 61):
        n = sum(1 for k in a if (k + off) in b and a[k][0] == b[k + off][0])
        if n > best_n:
            best, best_n = off, n
    return best, best_n


def kabsch(mobile, target):
    cm, ct = mobile.mean(0), target.mean(0)
    P, Q = mobile - cm, target - ct
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, ct - R @ cm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="models")
    ap.add_argument("--crystals", default="crystals")
    ap.add_argument("--results", default="results")
    a = ap.parse_args()

    # ---------- plot 1: per-residue deviation, model vs its own crystal -----
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    any_done = False
    for v in VARIANTS:
        mp = os.path.join(a.models, v + ".pdb")
        cp = os.path.join(a.crystals, v + ".pdb")
        if not (os.path.exists(mp) and os.path.exists(cp)):
            print(f"  {v}: missing a file, skipped")
            continue
        m, c = read_ca(mp), read_ca(cp)
        off, n_match = best_offset(m, c)
        keys = sorted(k for k in m if (k + off) in c and m[k][0] == c[k + off][0])
        if len(keys) < 50:
            print(f"  {v}: only {len(keys)} matching residues, skipped")
            continue
        P = np.array([m[k][1] for k in keys])
        Q = np.array([c[k + off][1] for k in keys])
        R, t = kabsch(P, Q)
        dev = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
        ax.plot(keys, dev, lw=1.1, label=f"{v}  (RMSD {np.sqrt((dev**2).mean()):.2f} A)")
        print(f"  {v:<6} {len(keys)} residues matched (offset {off:+d})   "
              f"RMSD {np.sqrt((dev**2).mean()):.2f} A   max {dev.max():.2f} A")
        any_done = True

    if any_done:
        for r in SITE:
            ax.axvline(r, color="0.78", ls=":", lw=0.8, zorder=0)
        ax.set_xlabel("Residue number")
        ax.set_ylabel(r"C$\alpha$ deviation from crystal ($\AA$)")
        ax.set_title("Where each AlphaFold2 model differs from its own crystal structure\n"
                     "(dotted lines mark active-site residues)", fontsize=10)
        ax.legend(fontsize=8, ncol=3)
        fig.tight_layout()
        fig.savefig("dev_per_residue.png", dpi=200)
        print("\nwrote dev_per_residue.png")
    plt.close(fig)

    # ---------- plot 2: active-site error, straight from results.json -------
    rj = os.path.join(a.results, "results.json")
    if not os.path.exists(rj):
        print(f"\n{rj} not found, skipping the active-site error plot")
        return
    S = json.load(open(rj))["structures"]

    present = [v for v in VARIANTS
               if f"AF2:{v}" in S and f"XTAL:{v}" in S]
    width = 0.8 / max(1, len(present))
    x = np.arange(len(SITE))
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for i, v in enumerate(present):
        errs = []
        for rn in SITE:
            af = S[f"AF2:{v}"]["site"].get(str(rn))
            xt = S[f"XTAL:{v}"]["site"].get(str(rn))
            errs.append(af["dist"] - xt["dist"] if af and xt else np.nan)
        ax.bar(x + i * width - 0.4 + width / 2, errs, width, label=v)
    ax.axhline(0, color="k", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{rn}" for rn in SITE])
    ax.set_xlabel("Active-site residue")
    ax.set_ylabel(r"model $-$ crystal distance to ligand ($\AA$)")
    ax.set_title("Error in where AlphaFold2 places each active-site residue\n"
                 "(negative means the model puts it too close to the ligand)",
                 fontsize=10)
    ax.legend(fontsize=8, ncol=5)
    fig.tight_layout()
    fig.savefig("site_error.png", dpi=200)
    plt.close(fig)
    print("wrote site_error.png")

    print("\nresidue 95, model minus crystal:")
    for v in present:
        af = S[f"AF2:{v}"]["site"].get("95")
        xt = S[f"XTAL:{v}"]["site"].get("95")
        if af and xt:
            print(f"  {v:<6} {af['aa']}95  model {af['dist']:.2f}  "
                  f"crystal {xt['dist']:.2f}  error {af['dist']-xt['dist']:+.2f} A")


if __name__ == "__main__":
    main()