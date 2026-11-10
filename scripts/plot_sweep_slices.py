#!/usr/bin/env python3
"""Central slice of selected reconstructions of the lambda sweep, side by side.

Rows: m_x, m_y, m_z. Columns: the given stems, in order. Every recon is divided by its
mean |m| inside the --mask support (TV shrinks the magnitude), and all panels share one
symmetric colour scale. Panels have no text (figstyle); the column order is printed.

Usage: plot_sweep_slices.py out_prefix stem [stem ...] [--mask mask.tiff] [--roi y0 x0 size]
Writes out_prefix_full.png (whole slice) and out_prefix_roi.png (ROI, default 128 px at the
centre of the mask).
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import numpy as np
import tifffile

from figstyle import image_grid


def load(stem):
    return np.stack([tifffile.imread(f"{stem}{i}.tiff").astype(np.float64) for i in range(3)])


def centre_crop(a, shape):
    oy, ox = (a.shape[0] - shape[0]) // 2, (a.shape[1] - shape[1]) // 2
    return a[oy : oy + shape[0], ox : ox + shape[1]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("prefix")
    ap.add_argument("stems", nargs="+")
    ap.add_argument("--mask", default="recons/emma/665/mask_tv.tiff")
    ap.add_argument("--roi", nargs=3, type=int, metavar=("Y0", "X0", "SIZE"))
    ap.add_argument("--vmax", type=float, default=1.3)
    ap.add_argument("--scale", type=int, default=3, help="output pixels per image pixel (full)")
    args = ap.parse_args()

    slices = []
    for s in args.stems:
        v = load(s)
        k = v.shape[1] // 2
        m2 = centre_crop(tifffile.imread(args.mask) > 0.5, v.shape[2:])
        mean_mag = np.linalg.norm(v[:, k], axis=0)[m2].mean()
        slices.append(v[:, k] / mean_mag)
        print(f"{s}: slice {k} of {v.shape[1]}, mean |m| {mean_mag:.4g}")
    # recons differ by a pixel or two in extent: centre-crop all to the common shape
    ny, nx = (min(s.shape[d] for s in slices) for d in (1, 2))
    slices = [s[:, (s.shape[1] - ny) // 2 :][:, :ny, (s.shape[2] - nx) // 2 :][:, :, :nx]
              for s in slices]
    m2 = centre_crop(tifffile.imread(args.mask) > 0.5, (ny, nx))

    if args.roi:
        y0, x0, size = args.roi
    else:
        ys, xs = np.nonzero(m2)
        size = 128
        y0, x0 = int(ys.mean()) - size // 2, int(xs.mean()) - size // 2

    for tag, (ya, xa, h, w, sc) in {
        "full": (0, 0, ny, nx, args.scale),
        "roi": (y0, x0, size, size, max(1, 1200 // (size * len(slices)) + 1) * 2),
    }.items():
        fig, axs = image_grid(3, len(slices), w * sc, h * sc)
        for j, sl in enumerate(slices):
            for i in range(3):
                axs[i, j].imshow(sl[i, ya : ya + h, xa : xa + w], cmap="RdBu_r",
                                 vmin=-args.vmax, vmax=args.vmax, interpolation="nearest")
        fig.savefig(f"{args.prefix}_{tag}.png", dpi=100)
        print(f"wrote {args.prefix}_{tag}.png ({tag}: y0 {ya}, x0 {xa}, {h}x{w})")


if __name__ == "__main__":
    main()
