#!/usr/bin/env python3
"""
eizs_pairwise.py -- Does AlphaFold2 reproduce the active-site remodelling that
                    a point mutation actually causes in epi-isozizaene synthase?

DESIGN


    delta(experiment) = pocket(mutant crystal) - pocket(WT crystal)
    delta(AlphaFold2) = pocket(mutant model)   - pocket(WT model)


Every structure is superposed by matched CA atoms
onto the WT crystal (3KB9), and ONE ligand set (3KB9's Mg3 . PPi . BTAC) is
used for all of them. If each crystal kept its own native
ligands while the models got transplanted ones, delta(experiment) would contain
ligand movement that delta(AlphaFold2) could not contain, and the two would not
be the same quantity. Holding the ligand fixed makes every volume difference a
difference in PROTEIN atoms.

This is not docking and predicts no binding pose. The ligand is a fixed
reference point defining where the active site is, identically everywhere. The
script reports, for every structure, the superposition RMSD and the minimum
ligand-to-protein distance, so the assumption is checked rather than asserted:
a clash means the transplant failed there and the number must not be used.

Usage

    python eizs_pairwise.py --models models/ --crystals crystals/ --out results/

Both directories hold one PDB per variant with MATCHING names, e.g.
    models/WT.pdb    models/F95H.pdb    models/F95M.pdb ...
    crystals/WT.pdb  crystals/F95H.pdb  crystals/F95M.pdb ...
"""

import argparse
import glob
import json
import os
import sys

import numpy as np

# constants

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

# Zamyatnin mean residue volumes (A^3) -- the null model: what the pocket would
# gain or lose from the substitution alone, with nothing else moving.
RES_VOL = {
    "G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7,
    "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2,
    "M": 162.9, "I": 166.7, "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9,
    "Y": 193.6, "W": 227.8,
}

ACTIVE_SITE = [69, 95, 96, 198, 203, 325]

SOLVENT = {"HOH", "WAT", "DOD", "SOL"}
EXCLUDE_HET = {
    "GOL", "EDO", "PEG", "PGE", "SO4", "PO4", "CL", "NA", "K", "ACT", "ACY",
    "DMS", "MPD", "TRS", "IMD", "FMT", "CIT", "EPE", "MES", "NO3", "ZN", "CA",
    "IOD", "BR", "CAC", "P6G", "1PE", "2PE", "PG4", "BME", "TLA",
}

VDW = 1.70
PROBE = 1.40
CONTACT = 4.5
CLASH = 2.2   # below this, a transplanted ligand atom overlaps protein


# parsing

class Structure:
    def __init__(self, path, label=None):
        self.path = path
        self.label = label or os.path.splitext(os.path.basename(path))[0]

        xyz, residx = [], []
        residues, order, het = {}, [], {}

        with open(path) as fh:
            for line in fh:
                rec = line[:6]
                if rec not in ("ATOM  ", "HETATM"):
                    continue
                if line[16] not in (" ", "A"):
                    continue
                resname = line[17:20].strip()
                atname = line[12:16].strip()
                element = line[76:78].strip().upper() or (
                    atname[0] if atname and not atname[0].isdigit() else "C")
                if element in ("H", "D"):
                    continue
                try:
                    p = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue

                is_aa = (rec == "ATOM  " and resname in AA3TO1) or resname == "MSE"
                if is_aa:
                    try:
                        seqnum = int(line[22:26])
                    except ValueError:
                        continue
                    key = (line[21], seqnum, line[26])
                    if key not in residues:
                        residues[key] = {"resseq": seqnum, "resname": resname,
                                         "aa": AA3TO1.get(resname, "X"),
                                         "ca": None, "bfac": 0.0, "atoms": []}
                        order.append(key)
                    r = residues[key]
                    r["atoms"].append(p)
                    if atname == "CA":
                        r["ca"] = p
                        try:
                            r["bfac"] = float(line[60:66])
                        except ValueError:
                            pass
                    xyz.append(p)
                    residx.append(len(order) - 1)
                elif rec == "HETATM" and resname not in SOLVENT \
                        and resname not in EXCLUDE_HET:
                    het.setdefault(resname, []).append(p)

        self.residues = [residues[k] for k in order]
        self.xyz = np.asarray(xyz, float).reshape(-1, 3)
        self.residx = np.asarray(residx, int)
        self.het = {k: np.asarray(v, float).reshape(-1, 3) for k, v in het.items()}
        self.seq = "".join(r["aa"] for r in self.residues)

    def ca(self, idx):
        return np.asarray([self.residues[i]["ca"] for i in idx], float)

    def by_resseq(self):
        return {r["resseq"]: i for i, r in enumerate(self.residues)}

    def __repr__(self):
        return (f"{self.label}: {len(self.residues)} res, {len(self.xyz)} atoms"
                + (f", het={sorted(self.het)}" if self.het else ""))


# alignment & superposition

def align(a, b):
    """Residue index pairs (a_i, b_i) by global sequence alignment.

    Used instead of matching residue numbers directly: the crystallised
    construct carries an expression tag and has disordered gaps, so author
    numbering cannot be trusted to line up with the model's 1..N numbering.
    """
    try:
        from Bio.Align import PairwiseAligner
    except ImportError:
        sys.exit("biopython required:  pip install biopython")
    al = PairwiseAligner()
    al.mode = "global"
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.match_score, al.mismatch_score = 2, -1
    best = al.align(a.seq, b.seq)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(best.aligned[0], best.aligned[1]):
        pairs += [(a0 + k, b0 + k) for k in range(a1 - a0)]
    return [(i, j) for i, j in pairs
            if a.residues[i]["ca"] is not None and b.residues[j]["ca"] is not None]


def kabsch(mobile, target):
    mobile, target = np.asarray(mobile, float), np.asarray(target, float)
    cm, ct = mobile.mean(0), target.mean(0)
    P, Q = mobile - cm, target - ct
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    fit = (R @ P.T).T
    rmsd = float(np.sqrt(np.mean(np.sum((fit - Q) ** 2, 1))))
    return R, ct - R @ cm, rmsd


def move(coords, R, t):
    c = np.asarray(coords, float)
    return c if len(c) == 0 else (R @ c.T).T + t


#  measurements

def pocket_volume(prot, centre, radius=9.0, spacing=0.6):
    """Connected free volume within `radius` of the pocket centre.

    Grid points further than vdW+probe from every protein heavy atom, keeping
    only the component containing the centre. The sphere cap makes the number
    directly comparable between structures and stops it collapsing to zero if
    a structure leaves the pocket mouth open.
    """
    from scipy.spatial import cKDTree
    from scipy import ndimage
    n = int(round(2 * radius / spacing)) + 1
    ax = np.linspace(-radius, radius, n)
    gx, gy, gz = np.meshgrid(ax, ax, ax, indexing="ij")
    local = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)
    inside = np.linalg.norm(local, axis=1) <= radius
    d, _ = cKDTree(prot).query(local + np.asarray(centre, float), k=1)
    mask = ((d > VDW + PROBE) & inside).reshape(n, n, n)
    if not mask.any():
        return 0.0
    lab, nlab = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    if nlab == 0:
        return 0.0
    c = n // 2
    seed = int(lab[c, c, c])
    if seed == 0:
        idx = np.argwhere(mask)
        i, j, k = idx[int(np.argmin(np.linalg.norm(idx - c, axis=1)))]
        seed = int(lab[i, j, k])
    return float((lab == seed).sum()) * spacing ** 3


def ligand_geometry(s, lig):
    """Min ligand-protein distance (clash check) + per-residue contacts."""
    from scipy.spatial import cKDTree
    d, _ = cKDTree(lig).query(s.xyz, k=1)
    per = {}
    for i, dist in enumerate(d):
        r = int(s.residx[i])
        rec = per.setdefault(r, {"min": 9e9, "n": 0})
        rec["min"] = min(rec["min"], float(dist))
        if dist <= CONTACT:
            rec["n"] += 1
    contacts = sorted(
        ({"resseq": s.residues[r]["resseq"], "aa": s.residues[r]["aa"],
          "n": v["n"], "min": round(v["min"], 2)} for r, v in per.items() if v["n"]),
        key=lambda x: (-x["n"], x["min"]))
    return float(d.min()), contacts


# main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True, help="AlphaFold2 models, <VARIANT>.pdb")
    ap.add_argument("--crystals", required=True, help="crystal structures, <VARIANT>.pdb")
    ap.add_argument("--out", default="results")
    ap.add_argument("--wt", default="WT")
    ap.add_argument("--site", default=",".join(map(str, ACTIVE_SITE)))
    ap.add_argument("--keep-het", default="")
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args()

    for c in (x.strip().upper() for x in a.keep_het.split(",") if x.strip()):
        EXCLUDE_HET.discard(c)
    os.makedirs(a.out, exist_ok=True)
    site = [int(x) for x in a.site.split(",") if x.strip()]

    #  reference frame: the WT crystal, and its ligands 
    ref_path = os.path.join(a.crystals, a.wt + ".pdb")
    if not os.path.exists(ref_path):
        sys.exit(f"WT crystal not found: {ref_path}")
    ref = Structure(ref_path, "XTAL:" + a.wt)
    if not ref.het:
        sys.exit("no ligands in the WT crystal -- did you download the closed "
                 "complex (3KB9)? try --keep-het BTM,POP,MG")
    lig = np.vstack(list(ref.het.values()))
    centre = lig.mean(0)
    print(f"reference  {ref}")
    print(f"ligand set " + ", ".join(f"{k}({len(v)})" for k, v in ref.het.items())
          + f"  centroid {np.round(centre, 1)}\n")

    # load everything, superpose onto the reference
    records = {}
    for kind, folder in (("AF2", a.models), ("XTAL", a.crystals)):
        for path in sorted(glob.glob(os.path.join(folder, "*.pdb"))):
            name = os.path.splitext(os.path.basename(path))[0]
            label = f"{kind}:{name}"
            s = Structure(path, label)
            pairs = align(s, ref)
            if len(pairs) < 50:
                print(f"  !! {label}: only {len(pairs)} aligned residues, skipped")
                continue
            R, t, rmsd = kabsch(s.ca([i for i, _ in pairs]),
                                ref.ca([j for _, j in pairs]))
            s.xyz = move(s.xyz, R, t)
            for r in s.residues:
                r["atoms"] = move(r["atoms"], R, t)
                if r["ca"] is not None:
                    r["ca"] = tuple(move([r["ca"]], R, t)[0])

            vol = pocket_volume(s.xyz, centre)
            mind, contacts = ligand_geometry(s, lig)
            idx = s.by_resseq()
            per_site = {}
            for rn in site:
                if rn in idx:
                    at = np.asarray(s.residues[idx[rn]]["atoms"], float)
                    per_site[rn] = {
                        "aa": s.residues[idx[rn]]["aa"],
                        "dist": round(float(np.min(np.linalg.norm(
                            at[:, None, :] - lig[None, :, :], axis=-1))), 2)}

            records[label] = {"kind": kind, "variant": name, "n_res": len(s.residues),
                              "rmsd_to_ref": round(rmsd, 2), "pocket": round(vol, 1),
                              "min_lig_dist": round(mind, 2),
                              "n_contacts": len(contacts), "contacts": contacts[:25],
                              "site": per_site}
            flag = "  <-- CLASH" if mind < CLASH else ""
            print(f"  {label:<14} RMSD {rmsd:5.2f} A   pocket {vol:6.0f} A^3   "
                  f"min-lig {mind:4.2f} A   contacts {len(contacts):2d}{flag}")

    #  the pairwise comparison 
    wt_af2, wt_xt = f"AF2:{a.wt}", f"XTAL:{a.wt}"
    for need in (wt_af2, wt_xt):
        if need not in records:
            sys.exit(f"missing {need} -- both WT model and WT crystal are required")

    variants = sorted({r["variant"] for r in records.values()} - {a.wt})
    pairs_out = []
    for v in variants:
        m, x = records.get(f"AF2:{v}"), records.get(f"XTAL:{v}")
        if not (m and x):
            continue
        d_exp = x["pocket"] - records[wt_xt]["pocket"]
        d_af2 = m["pocket"] - records[wt_af2]["pocket"]
        dv = None
        if len(v) >= 3 and v[0].isalpha() and v[-1].isalpha():
            try:
                f, to = v[0].upper(), v[-1].upper()
                int(v[1:-1])
                if f in RES_VOL and to in RES_VOL:
                    dv = round(RES_VOL[to] - RES_VOL[f], 1)
            except ValueError:
                pass
        pairs_out.append({"variant": v, "d_exp": round(d_exp, 1),
                          "d_af2": round(d_af2, 1),
                          "residual": round(d_af2 - d_exp, 1),
                          "dV_sidechain": dv,
                          "same_sign": bool(d_exp * d_af2 > 0)})

    print("\n  pairwise:")
    for p in pairs_out:
        print(f"    {p['variant']:<8} d(exp) {p['d_exp']:+7.1f}   "
              f"d(AF2) {p['d_af2']:+7.1f}   residual {p['residual']:+7.1f} A^3"
              + ("   same sign" if p["same_sign"] else "   OPPOSITE SIGN"))

    summary = {"reference": os.path.basename(ref_path), "site": site,
               "ligands": {k: len(v) for k, v in ref.het.items()},
               "structures": records, "pairs": pairs_out}
    with open(os.path.join(a.out, "results.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    write_md(summary, a, os.path.join(a.out, "results.md"))
    if not a.no_figures:
        try:
            figures(summary, a)
        except Exception as e:
            print(f"  [warn] figures failed: {e}")
    print(f"\nwrote {a.out}/results.md and results.json")


def write_md(S, a, path):
    R, L = S["structures"], []
    order = sorted(R, key=lambda k: (R[k]["variant"] != a.wt, R[k]["variant"], R[k]["kind"]))

    L.append("# Does AlphaFold2 reproduce the real active-site change?\n")
    L.append(f"Reference frame and ligand source: **{S['reference']}** "
             + ", ".join(f"{k} ({v} atoms)" for k, v in S["ligands"].items())
             + ". Every structure is superposed onto it and shares this one "
               "ligand set, so all volume differences are differences in protein atoms.\n")

    L.append("\n## Table 1 - Every structure\n")
    L.append("| Structure | Residues | CA RMSD to ref (A) | Pocket volume (A^3) "
             "| Min ligand-protein distance (A) | Contact residues |")
    L.append("|---|---|---|---|---|---|")
    for k in order:
        r = R[k]
        warn = " **CLASH**" if r["min_lig_dist"] < CLASH else ""
        L.append(f"| {k} | {r['n_res']} | {r['rmsd_to_ref']:.2f} | {r['pocket']:.0f} "
                 f"| {r['min_lig_dist']:.2f}{warn} | {r['n_contacts']} |")
    L.append(f"\nAny value below {CLASH} A in the fifth column means the transplanted "
             "ligand overlaps protein in that structure and its pocket volume must "
             "not be used.\n")

    L.append("\n## Table 2 - THE RESULT: predicted vs real change\n")
    L.append("| Mutant | dV side chain (A^3) | d Pocket, experiment (A^3) "
             "| d Pocket, AlphaFold2 (A^3) | Residual (A^3) | Same sign? |")
    L.append("|---|---|---|---|---|---|")
    for p in S["pairs"]:
        dv = f"{p['dV_sidechain']:+.1f}" if p["dV_sidechain"] is not None else "-"
        L.append(f"| {p['variant']} | {dv} | {p['d_exp']:+.1f} | {p['d_af2']:+.1f} "
                 f"| {p['residual']:+.1f} | {'yes' if p['same_sign'] else '**no**'} |")
    L.append("\n*d Pocket, experiment* = mutant crystal minus WT crystal: what the "
             "mutation really does. *d Pocket, AlphaFold2* = mutant model minus WT "
             "model: what AlphaFold2 predicts from sequence alone. *Residual* is the "
             "gap between them. *dV side chain* is the null model -- the volume the "
             "substitution adds or removes with nothing else moving.\n")

    L.append("\n## Table 3 - Minimum residue-ligand distance (A)\n")
    L.append("| Structure | " + " | ".join(str(s) for s in S["site"]) + " |")
    L.append("|---" * (len(S["site"]) + 1) + "|")
    for k in order:
        cells = []
        for rn in S["site"]:
            e = R[k]["site"].get(str(rn), R[k]["site"].get(rn))
            cells.append(f"{e['aa']}{rn}: {e['dist']:.2f}" if e else "n/a")
        L.append(f"| {k} | " + " | ".join(cells) + " |")
    L.append("\nCheck the XTAL:WT row shows F95, F96, F198, W203. If it does not, the "
             "crystal numbering does not match and nothing above is valid.\n")

    L.append(f"\n## Table 4 - Residues contacting the ligands, {a.wt} crystal\n")
    L.append("| Residue | Contacts | Min distance (A) |")
    L.append("|---|---|---|")
    for c in R[f"XTAL:{a.wt}"]["contacts"][:20]:
        L.append(f"| {c['aa']}{c['resseq']} | {c['n']} | {c['min']:.2f} |")

    open(path, "w").write("\n".join(L) + "\n")


def figures(S, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    P = S["pairs"]
    if P:
        fig, ax = plt.subplots(figsize=(5.2, 5))
        xs = [p["d_exp"] for p in P]
        ys = [p["d_af2"] for p in P]
        lim = max(10, max(abs(v) for v in xs + ys) * 1.25)
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="perfect agreement")
        ax.axhline(0, color="0.7", lw=0.6)
        ax.axvline(0, color="0.7", lw=0.6)
        ax.scatter(xs, ys, s=70, zorder=3, color="#3b76af")
        for p in P:
            ax.annotate(p["variant"], (p["d_exp"], p["d_af2"]),
                        textcoords="offset points", xytext=(7, 5), fontsize=9)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel(r"$\Delta$ pocket volume, experiment ($\AA^3$)")
        ax.set_ylabel(r"$\Delta$ pocket volume, AlphaFold2 ($\AA^3$)")
        ax.set_title("Predicted vs real active-site change")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig1_predicted_vs_real.png"), dpi=200)
        plt.close(fig)

    R = S["structures"]
    variants = [a.wt] + sorted({r["variant"] for r in R.values()} - {a.wt})
    xt = [R.get(f"XTAL:{v}", {}).get("pocket", np.nan) for v in variants]
    af = [R.get(f"AF2:{v}", {}).get("pocket", np.nan) for v in variants]
    x = np.arange(len(variants))
    fig, ax = plt.subplots(figsize=(7.5, 4))
    ax.bar(x - 0.2, xt, 0.4, label="crystal", color="#444")
    ax.bar(x + 0.2, af, 0.4, label="AlphaFold2", color="#3b76af")
    ax.set_xticks(x)
    ax.set_xticklabels(variants)
    ax.set_ylabel(r"Pocket free volume ($\AA^3$)")
    ax.set_title("Active-site pocket volume: experiment vs prediction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig2_pocket_volumes.png"), dpi=200)
    plt.close(fig)

    site = S["site"]
    labels = sorted(R, key=lambda k: (R[k]["variant"] != a.wt, R[k]["variant"], R[k]["kind"]))
    M = np.full((len(labels), len(site)), np.nan)
    for i, k in enumerate(labels):
        for j, rn in enumerate(site):
            e = R[k]["site"].get(str(rn), R[k]["site"].get(rn))
            if e:
                M[i, j] = e["dist"]
    fig, ax = plt.subplots(figsize=(6.5, 0.45 * len(labels) + 1.8))
    im = ax.imshow(M, cmap="viridis_r", aspect="auto")
    ax.set_xticks(range(len(site)))
    ax.set_xticklabels([str(s) for s in site])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(len(labels)):
        for j in range(len(site)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center",
                        color="w", fontsize=7)
    ax.set_xlabel("Active-site residue")
    ax.set_title(r"Minimum residue$-$ligand distance ($\AA$)")
    fig.colorbar(im, ax=ax, label=r"$\AA$")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig3_contacts.png"), dpi=200)
    plt.close(fig)
    print("  figures: fig1_predicted_vs_real.png, fig2_pocket_volumes.png, fig3_contacts.png")


if __name__ == "__main__":
    main()