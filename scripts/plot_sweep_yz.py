#!/usr/bin/env python3
"""y-z cross-section (fixed x) of selected reconstructions, side by side.

Rows: the given stems, in order. Columns: m_x, m_y, m_z. Horizontal axis y (--length pixels
from --y0), vertical axis z (all slices, centred on a common extent: a recon with fewer
slices is blank outside its own thickness). Same normalisation and colour scale as
plot_sweep_slices.py (division by the mean |m| of the central slice inside --mask).
Panels have no text (figstyle).

Usage: plot_sweep_yz.py out.png stem [stem ...] [--mask mask.tiff] [--x X] [--y0 Y0]
                        [--length 200] [--scale 4] [--zstretch 3]
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import numpy as np
import tifffile

from figstyle import image_grid
from plot_sweep_slices import centre_crop, load


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out")
    ap.add_argument("stems", nargs="+")
    ap.add_argument("--mask", default="recons/emma/665/mask_tv.tiff")
    ap.add_argument("--x", type=int, help="column of the section (default: centre of the support)")
    ap.add_argument("--y0", type=int, help="first row (default: centred on the support)")
    ap.add_argument("--length", type=int, default=200)
    ap.add_argument("--scale", type=int, default=4, help="output pixels per image pixel in y")
    ap.add_argument("--zstretch", type=int, default=3, help="extra output pixels per slice factor in z")
    ap.add_argument("--vmax", type=float, default=1.3)
    args = ap.parse_args()

    vols = []
    for s in args.stems:
        v = load(s)
        k = v.shape[1] // 2
        m2 = centre_crop(tifffile.imread(args.mask) > 0.5, v.shape[2:])
        vols.append(v / np.linalg.norm(v[:, k], axis=0)[m2].mean())
    ny, nx = (min(v.shape[d] for v in vols) for d in (2, 3))
    nzmax = max(v.shape[1] for v in vols)
    vols = [v[:, :, (v.shape[2] - ny) // 2 :][:, :, :ny, (v.shape[3] - nx) // 2 :][..., :nx]
            for v in vols]
    m2 = centre_crop(tifffile.imread(args.mask) > 0.5, (ny, nx))
    ys, xs = np.nonzero(m2)
    x = args.x if args.x is not None else int(xs.mean())
    y0 = args.y0 if args.y0 is not None else int(ys.mean()) - args.length // 2
    sl = slice(y0, y0 + args.length)
    print(f"section x = {x}, y = {y0}..{y0 + args.length - 1}, {nzmax} slices; stems: "
          + ", ".join(f"{s.split('/')[-2]}({v.shape[1]})" for s, v in zip(args.stems, vols)))

    fig, axs = image_grid(len(vols), 3, args.length * args.scale,
                          nzmax * args.scale * args.zstretch)
    for j, v in enumerate(vols):
        pad = (nzmax - v.shape[1]) // 2
        for i in range(3):
            img = np.full((nzmax, args.length), np.nan)
            img[pad : pad + v.shape[1]] = v[i, :, sl, x]
            axs[j, i].imshow(img, cmap="RdBu_r", vmin=-args.vmax, vmax=args.vmax,
                             interpolation="nearest", aspect="auto")
    fig.savefig(args.out, dpi=100)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
