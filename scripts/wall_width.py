#!/usr/bin/env python3
"""Domain-wall width of vector reconstructions.

Works on the direction field mh = m/|m|, slice by slice. Two estimators:

  delta_fit  Wall centre lines (skeleton of the thresholded gradient energy) are found; at
             every point mh is sampled along the wall normal and p(s) = mh(s).u, with u the
             unit vector from mh(-L) to mh(+L), is fitted by a*tanh((s-s0)/delta)+c.
             Accepted if the wall turns mh by more than 120 deg, R^2 > 0.9 and the energy
             profile has a single peak. Reported: median delta, IQR, number of walls.
  delta_g    2*TV/(pi*D) with TV = sum|grad mh|, D = sum|grad mh|^2 over a band around the
             walls. Exact for a tanh wall (|theta'| = sech(s/delta)/delta); no fit, but
             noise inflates D and biases it low.

Also the noise inside the domains (rms angle of mh to its local mean direction, away from
walls). With several reconstructions, the wall points of --ref (default: the first) are
fitted in all of them (paired comparison) and the "common" columns use the walls accepted
in every reconstruction.

Reconstructions are <stem>{0,1,2}.tiff (volume frame, centre-cropped mask), or --ovf files.

Usage: wall_width.py stem [stem ...] [--mask mask.tiff] [--thickness N] [--ref stem]
                     [--half-length px] [--pixel-size nm] [--png out.png] [--json out.json]
       wall_width.py --ovf truth.ovf [stem ...] --mask-none
"""

import argparse
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
from skimage.filters import threshold_otsu
from skimage.morphology import skeletonize

L = 15  # half length of the profiles [px]
R2_MIN = 0.9


def load_tiff(stem):
    return np.stack([tifffile.imread(f"{stem}{i}.tiff").astype(np.float64) for i in range(3)])


def load_ovf(path):
    with open(path) as f:
        text = f.read()
    head, _, rest = text.partition("# Begin: Data Text\n")
    hdr = {}
    for line in head.splitlines():
        if line.startswith("#") and ":" in line:
            k, v = line[1:].split(":", 1)
            hdr[k.strip()] = v.strip()
    nx, ny, nz = (int(hdr[k]) for k in ("xnodes", "ynodes", "znodes"))
    data = np.array(rest.split("# End: Data Text")[0].split(), dtype=np.float64)
    return np.moveaxis(data.reshape(nz, ny, nx, 3), -1, 0)


def centre_crop(a, shape):
    oy, ox = (a.shape[0] - shape[0]) // 2, (a.shape[1] - shape[1]) // 2
    return a[oy : oy + shape[0], ox : ox + shape[1]]


def unit(m):
    return m / np.maximum(np.linalg.norm(m, axis=0), 1e-12)


def grads(mh, sigma=1.0):
    """d(mh_i)/d(y,x) of one slice mh (3, ny, nx) -> (3, 2, ny, nx), Gaussian-smoothed if sigma."""
    s = np.stack([ndi.gaussian_filter(c, sigma) for c in mh]) if sigma else mh
    return np.stack([np.stack(np.gradient(c)) for c in s])


def detect(mh, mask):
    """Wall points of one slice: list of (y, x, ny, nx) and the binary wall band."""
    g = grads(mh)
    e = (g**2).sum(axis=(0, 1))
    e[~mask] = 0
    if not mask.any() or e[mask].max() <= 0:
        return [], np.zeros_like(mask)
    wall = e > threshold_otsu(e[mask])
    wall &= mask
    skel = skeletonize(wall)
    # structure tensor of the gradients: the leading eigenvector is the wall normal
    jyy = ndi.gaussian_filter((g[:, 0] * g[:, 0]).sum(0), 2.0)
    jxx = ndi.gaussian_filter((g[:, 1] * g[:, 1]).sum(0), 2.0)
    jxy = ndi.gaussian_filter((g[:, 0] * g[:, 1]).sum(0), 2.0)
    ang = 0.5 * np.arctan2(2 * jxy, jxx - jyy)  # direction of the largest eigenvalue (x, y)
    ys, xs = np.nonzero(skel)
    pts = [(y, x, np.sin(ang[y, x]), np.cos(ang[y, x])) for y, x in zip(ys, xs)]
    return pts, wall


def tanh(s, a, s0, d, c):
    return a * np.tanh((s - s0) / d) + c


def fit_point(m, pt, normalize=True):
    """Fit the wall at pt = (y, x, ny, nx) in one slice of the field m (3, ny, nx).

    With normalize the profile is the direction field; otherwise the raw vector (scaled to
    mean |m| = 1 by the caller), which avoids the sharpening caused by |m| dipping in the wall.
    -> (delta, R2) or None."""
    y, x, ny, nx = pt
    s = np.arange(-L, L + 1, dtype=float)
    yy, xx = y + s * ny, x + s * nx
    prof = np.stack([ndi.map_coordinates(c, [yy, xx], order=1, mode="nearest") for c in m])
    prof_u = prof / np.maximum(np.linalg.norm(prof, axis=0), 1e-12)
    if not normalize:
        prof_n = prof
        prof = prof_u  # shape tests use the direction
    else:
        prof_n = prof_u
    a, b = prof_u[:, 0], prof_u[:, -1]
    if a @ b > np.cos(np.deg2rad(120.0)):
        return None
    u = b - a
    u /= np.linalg.norm(u)
    p = u @ prof_n
    # single energy peak along the profile
    de = (np.diff(prof, axis=1) ** 2).sum(0)
    peaks = (de[1:-1] > de[:-2]) & (de[1:-1] >= de[2:]) & (de[1:-1] > 0.5 * de.max())
    if peaks.sum() > 1:
        return None
    try:
        popt, _ = curve_fit(tanh, s, p, p0=[1.0, 0.0, 3.0, 0.0], maxfev=400,
                            bounds=([0, -L, 0.2, -1], [2, L, 4 * L, 1]))
    except (RuntimeError, ValueError):
        return None
    res = p - tanh(s, *popt)
    r2 = 1 - (res @ res) / max(((p - p.mean()) ** 2).sum(), 1e-12)
    if r2 < R2_MIN or abs(popt[1]) > L / 2:
        return None
    return popt[2], r2


def band_metrics(mh, wall, mask):
    """delta_g and domain noise of one slice."""
    band = ndi.binary_dilation(wall, iterations=L) & mask
    g = grads(mh, sigma=0)  # unsmoothed: smoothing widens walls that are a few px thick
    gm = np.sqrt((g**2).sum(axis=(0, 1)))
    tv, d = gm[band].sum(), (gm[band] ** 2).sum()
    dom = mask & ~ndi.binary_dilation(wall, iterations=max(L, 4))
    noise = None
    if dom.sum() > 100:  # rms |grad mh| in the domains [deg/px]
        noise = float(np.rad2deg(np.sqrt(np.mean(gm[dom] ** 2))))
    return tv, d, noise


def analyse(mhs, names, masks, slices, ref, fields, normalize):
    """mhs: list of (3, nz, ny, nx) direction fields. Returns per-recon results."""
    ref_i = names.index(ref) if ref in names else 0
    res = [dict(name=n, deltas=[], tv=0.0, d=0.0, noise=[]) for n in names]
    fits = [[] for _ in names]  # per recon: list of (k, y, x, delta or None)
    for k in slices:
        pts, wall = detect(mhs[ref_i][:, k], masks[k])
        for i, mh in enumerate(mhs):
            if i != ref_i:
                _, w_i = detect(mh[:, k], masks[k])
            else:
                w_i = wall
            tv, d, noise = band_metrics(mh[:, k], w_i, masks[k])
            res[i]["tv"] += tv
            res[i]["d"] += d
            if noise is not None:
                res[i]["noise"].append(noise)
            for pt in pts:
                f = fit_point(fields[i][:, k], pt, normalize)
                fits[i].append((k, pt[0], pt[1], None if f is None else f[0]))
    return res, fits


def summarize(res, fits, px):
    common = np.array([all(f[j][3] is not None for f in fits) for j in range(len(fits[0]))])
    out = []
    for r, f in zip(res, fits):
        dl = np.array([t[3] if t[3] is not None else np.nan for t in f])
        ok = ~np.isnan(dl)
        row = dict(name=r["name"], n_points=len(f), n_ok=int(ok.sum()), n_common=int(common.sum()))
        for tag, sel in (("", ok), ("_common", common)):
            v = dl[sel]
            if v.size:
                q1, q2, q3 = np.percentile(v, [25, 50, 75])
                row[f"delta_fit{tag}"] = float(q2 * px)
                row[f"iqr{tag}"] = float((q3 - q1) * px)
        row["delta_g"] = float(2 * r["tv"] / (np.pi * r["d"]) * px) if r["d"] > 0 else None
        row["noise_deg"] = float(np.mean(r["noise"])) if r["noise"] else None
        out.append(row)
    return out


def main():
    global L
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stems", nargs="*")
    ap.add_argument("--ovf", help="ground-truth OVF, analysed as an extra entry named 'truth'")
    ap.add_argument("--mask", default="recons/emma/665/mask_tv.tiff")
    ap.add_argument("--mask-none", action="store_true", help="no support mask")
    ap.add_argument("--erode", type=int, default=6, help="mask erosion [px]")
    ap.add_argument("--thickness", type=int, default=5, help="N central slices")
    ap.add_argument("--ref", help="stem (or 'truth') whose walls are fitted everywhere")
    ap.add_argument("--pixel-size", type=float, default=1.0, help="nm per px (default: px)")
    ap.add_argument("--half-length", type=int, default=L,
                    help="half length of the profiles [px]; must stay below the domain width/2")
    ap.add_argument("--normalize", action="store_true",
                    help="fit the direction field m/|m| instead of the raw m (sharpens walls "
                    "where |m| dips)")
    ap.add_argument("--png")
    ap.add_argument("--json")
    args = ap.parse_args()
    L = args.half_length

    names, vols = [], []
    if args.ovf:
        names.append("truth")
        vols.append(load_ovf(args.ovf))
    for s in args.stems:
        names.append(s)
        vols.append(load_tiff(s))
    if not vols:
        sys.exit("nothing to analyse")
    # volumes of different extent (padding, mask_dims output): centre-crop to the common shape
    sh = [min(v.shape[d] for v in vols) for d in (1, 2, 3)]
    for i, v in enumerate(vols):
        o = [(n - m) // 2 for n, m in zip(v.shape[1:], sh)]
        vols[i] = v[:, o[0] : o[0] + sh[0], o[1] : o[1] + sh[1], o[2] : o[2] + sh[2]]
    if len({v.shape for v in vols}) != 1:
        sys.exit(f"shape mismatch: {[v.shape for v in vols]}")
    _, nz, ny, nx = vols[0].shape
    m2 = np.ones((ny, nx), bool) if args.mask_none else centre_crop(tifffile.imread(args.mask) > 0.5, (ny, nx))
    m2 = ndi.binary_erosion(m2, iterations=args.erode, border_value=0)
    masks = [m2] * nz
    z0 = (nz - args.thickness) // 2
    slices = range(max(z0, 0), min(z0 + args.thickness, nz))
    mhs = [unit(v) for v in vols]
    # raw fields scaled to mean |m| = 1 inside the support (the amplitude is a free parameter)
    fields = [v / np.linalg.norm(v, axis=0)[:, m2].mean() for v in vols]
    print(f"volume {nz}x{ny}x{nx}, slices {slices.start}..{slices.stop - 1}, "
          f"support fraction {m2.mean():.3f}", file=sys.stderr)

    ref = args.ref or names[0]
    res, fits = analyse(mhs, names, masks, slices, ref, mhs if args.normalize else fields, args.normalize)
    out = summarize(res, fits, args.pixel_size)
    unitname = "nm" if args.pixel_size != 1.0 else "px"
    print(f"walls detected in {ref}; widths in {unitname} (delta of tanh; Lilley width = pi*delta)")
    print(f"{'recon':50s} {'N':>6s} {'fit':>6s} {'IQR':>6s} {'common':>7s} {'IQRc':>6s} {'delta_g':>8s} {'noise':>6s}")
    for r in out:
        print(f"{r['name'][-50:]:50s} {r['n_ok']:6d}/{r['n_points']:<6d} "
              f"{r.get('delta_fit', float('nan')):6.2f} {r.get('iqr', float('nan')):6.2f} "
              f"{r.get('delta_fit_common', float('nan')):7.2f} {r.get('iqr_common', float('nan')):6.2f} "
              f"{r['delta_g'] if r['delta_g'] is not None else float('nan'):8.2f} "
              f"{r['noise_deg'] if r['noise_deg'] is not None else float('nan'):6.1f}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
    if args.png:
        k = (slices.start + slices.stop - 1) // 2
        fig, ax = plt.subplots(1, 2, figsize=(11, 5))
        ri = names.index(ref) if ref in names else 0
        ax[0].imshow(vols[ri][2, k], cmap="gray")
        pk = [(f[1], f[2], f[3]) for f in fits[ri] if f[0] == k]
        for y, x, d in pk:
            ax[0].plot(x, y, ".", ms=2, color="tab:green" if d is not None else "tab:red")
        ax[0].set_title(f"{ref[-40:]} m_z, slice {k}: accepted (green) / rejected (red)")
        for n, f in zip(names, fits):
            v = [t[3] for t in f if t[3] is not None]
            ax[1].hist(np.asarray(v) * args.pixel_size, bins=40, histtype="step", label=n[-30:])
        ax[1].set_xlabel(f"delta [{unitname}]")
        ax[1].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(args.png, dpi=130)


if __name__ == "__main__":
    main()
