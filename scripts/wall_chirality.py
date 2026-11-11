#!/usr/bin/env python3
"""Bloch / Neel character of domain walls from the integrated in-plane moment.

Wall points and normals come from the first entry (--ovf truth, or the first stem) via
wall_width.detect. At every point the raw in-plane moment (m_x, m_y) is summed along the
wall normal over s = -L..L (1 px steps):

  I_n = sum m_par . n      across the wall   (Neel-like)
  I_t = sum m_par . t      along the wall    (Bloch-like), t = z x n

n is oriented so that m_z increases along +n (each recon's own raw m_z). Blurring flattens
the in-plane peak but conserves this integral (the imaging is linear), unlike the peak or
the width. Units: Ms*px, Ms = mean |m| inside the support of the slice. For a tanh wall of
width delta with full in-plane amplitude Ms, |I| = pi*delta*Ms.

Reported per recon, over the walls accepted in all entries (m_z turns by > 1 Ms across the
profile in every entry): rms I_n, rms I_t, f_B = <I_t^2> / (<I_t^2> + <I_n^2>) (1 = Bloch,
0 = Neel, 0.5 = noise or mixed), and the noise floor, the robust rms (1.4826 MAD) of I over random profiles
(random in-plane direction) inside the domains away from the holes, same sampling. The corrected values are
sqrt(max(rms^2 - floor^2, 0)). Also the per-component Pearson r with the first entry (a
negative r would flip the chirality sign).

Usage: wall_chirality.py [--ovf truth.ovf] [stem ...] [--slices K ...] [--half-length 6]
                         [--mask-none | --mask mask.tiff] [--erode 6]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
from wall_width import centre_crop, detect, load_ovf, load_tiff, unit  # noqa: E402


def profiles(m, ys, xs, ny, nx, L):
    """Raw field m (3, Y, X) sampled along the normals: (3, N, 2L+1)."""
    s = np.arange(-L, L + 1, dtype=float)
    yy = ys[:, None] + s[None] * ny[:, None]
    xx = xs[:, None] + s[None] * nx[:, None]
    return np.stack([ndi.map_coordinates(c, [yy, xx], order=1, mode="nearest") for c in m])


def integrals(m, ys, xs, ny, nx, L):
    """-> (I_n, I_t, dmz, ok-by-oriented-sign) with n oriented along rising m_z."""
    p = profiles(m, ys, xs, ny, nx, L)
    k = max(2, L // 3)
    dmz = p[2][:, -k:].mean(1) - p[2][:, :k].mean(1)
    sgn = np.where(dmz >= 0, 1.0, -1.0)
    nxo, nyo = nx * sgn, ny * sgn  # oriented normal (x, y components)
    i_n = (p[0] * nxo[:, None] + p[1] * nyo[:, None]).sum(1)
    i_t = (-p[0] * nyo[:, None] + p[1] * nxo[:, None]).sum(1)
    return i_n, i_t, np.abs(dmz)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stems", nargs="*")
    ap.add_argument("--ovf")
    ap.add_argument("--mask", default="recons/emma/665/mask_tv.tiff")
    ap.add_argument("--mask-none", action="store_true")
    ap.add_argument("--erode", type=int, default=6)
    ap.add_argument("--slices", type=int, nargs="+", help="z indices (default: 5 central)")
    ap.add_argument("--half-length", type=int, default=6)
    ap.add_argument("--hole-margin", type=int, default=15,
                    help="px excluded around the truth's holes (|m| < 0.5)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json", help="write the table as a list of dicts")
    args = ap.parse_args()
    L = args.half_length

    names, vols = [], []
    if args.ovf:
        names.append("truth")
        vols.append(load_ovf(args.ovf))
    for s in args.stems:
        names.append(s)
        vols.append(load_tiff(s))
    sh = [min(v.shape[d] for v in vols) for d in (1, 2, 3)]
    for i, v in enumerate(vols):
        o = [(n - m) // 2 for n, m in zip(v.shape[1:], sh)]
        vols[i] = v[:, o[0] : o[0] + sh[0], o[1] : o[1] + sh[1], o[2] : o[2] + sh[2]]
    nz, ny_, nx_ = sh
    m2 = np.ones((ny_, nx_), bool) if args.mask_none else centre_crop(
        tifffile.imread(args.mask) > 0.5, (ny_, nx_))
    m2 = ndi.binary_erosion(m2, iterations=args.erode, border_value=0)
    if args.ovf:  # fiducial holes of the truth (|m| ~ 0) and the artefacts around them
        holes = np.linalg.norm(vols[0][:, nz // 2], axis=0) < 0.5
        m2 &= ~ndi.binary_dilation(holes, iterations=args.hole_margin)
    slices = args.slices or list(range((nz - 5) // 2, (nz - 5) // 2 + 5))
    print(f"volume {nz}x{ny_}x{nx_}, slices {slices}, L = {L}")

    rng = np.random.default_rng(args.seed)
    scale = [np.linalg.norm(v[:, nz // 2], axis=0)[m2].mean() for v in vols]
    # per entry: lists over all wall points of the reference
    res = [dict(i_n=[], i_t=[], dmz=[], noise_n=[], noise_t=[]) for _ in vols]
    for k in slices:
        pts, wall = detect(unit(vols[0][:, k]), m2)
        if not pts:
            continue
        ys, xs, nys, nxs = (np.array(c, float) for c in zip(*pts))
        dom = ndi.distance_transform_edt(~wall) > L + 2
        dom &= m2
        dy, dx = np.nonzero(dom)
        pick = rng.choice(len(dy), size=min(len(ys), len(dy)), replace=False)
        phi = rng.uniform(0, np.pi, size=len(pick))
        for i, v in enumerate(vols):
            m = v[:, k] / scale[i]
            a, b, d = integrals(m, ys, xs, nys, nxs, L)
            res[i]["i_n"].append(a)
            res[i]["i_t"].append(b)
            res[i]["dmz"].append(d)
            a, b, _ = integrals(m, dy[pick].astype(float), dx[pick].astype(float),
                                np.sin(phi), np.cos(phi), L)
            res[i]["noise_n"].append(a)
            res[i]["noise_t"].append(b)
    for r in res:
        for key in r:
            r[key] = np.concatenate(r[key])

    ok = np.all([r["dmz"] > 1.0 for r in res], axis=0)
    print(f"wall points: {len(ok)}, accepted in all entries: {int(ok.sum())}")
    rows = []
    print(f"{'entry':48s} {'N':>6s} {'rms In':>7s} {'rms It':>7s} {'floor':>6s} "
          f"{'In corr':>8s} {'It corr':>8s} {'f_B':>5s} {'f_B corr':>9s} {'r_xyz vs 1st'}")
    for i, (n, r) in enumerate(zip(names, res)):
        rn, rt = np.sqrt(np.mean(r["i_n"][ok] ** 2)), np.sqrt(np.mean(r["i_t"][ok] ** 2))
        mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))  # noqa: E731
        fl = 0.5 * (mad(r["noise_n"]) + mad(r["noise_t"]))  # robust: outliers excluded
        cn, ct = np.sqrt(max(rn**2 - fl**2, 0)), np.sqrt(max(rt**2 - fl**2, 0))
        fb = rt**2 / (rt**2 + rn**2)
        fbc = ct**2 / max(ct**2 + cn**2, 1e-12)
        k = nz // 2
        cc = [np.corrcoef(vols[0][c, k][m2], vols[i][c, k][m2])[0, 1] for c in range(3)]
        rows.append(dict(name=n, n=int(ok.sum()), rms_in=float(rn), rms_it=float(rt), floor=float(fl),
                         in_corr=float(cn), it_corr=float(ct), f_b=float(fb), f_b_corr=float(fbc),
                         r_x=float(cc[0]), r_y=float(cc[1]), r_z=float(cc[2])))
        print(f"{n[-48:]:48s} {int(ok.sum()):6d} {rn:7.2f} {rt:7.2f} {fl:6.2f} {cn:8.2f} "
              f"{ct:8.2f} {fb:5.2f} {fbc:9.2f}  " + " ".join(f"{c:+.2f}" for c in cc))
    if args.json:
        import json

        with open(args.json, "w") as f:
            json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main()
