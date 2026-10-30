#!/usr/bin/env python3
"""Add white Gaussian noise to projections, matched to the noise level of measured data.

The noise level of a projection stack is measured by the same high-pass estimator on the
central --window x --window pixels of every projection:

  h   = p - gaussian_filter(p, --smooth)          (removes the smooth signal)
  rho = 1.4826 * MAD(h) / std(p)                  (median over the projections)

The measured data (--like, several files are averaged) give rho_data; the clean simulated
stack gives rho_0 (the high-pass content of the signal itself: walls, edges). Noise of
relative std s = sigma / (std of the window) adds in quadrature, rho^2 (1 + s^2) = rho_0^2 + (0.946 s)^2 (0.946 is
the response of the estimator to white noise for --smooth 1.5), so

  s = sqrt((rho_data^2 - rho_0^2) / (0.946^2 - rho_data^2)) * --scale

which assumes the real signal is as sharp as the simulated one (otherwise the noise is
overestimated). rho_data is also printed per file. With --report-only nothing is written.

Usage: add_noise.py clean.tiff noisy.tiff --like data.tiff [data2.tiff ...] [--scale 1]
                    [--seed 0] [--report-only]
"""

import argparse

import numpy as np
import tifffile
from scipy import ndimage as ndi

WHITE_RESPONSE = 0.946  # std of (n - G*n) / std(n), 2-D Gaussian with sigma 1.5


def rho(stack, window, smooth):
    """-> (rho, median std of the projections over the central window)."""
    ny, nx = stack.shape[1:]
    y0, x0 = (ny - window) // 2, (nx - window) // 2
    out, sds = [], []
    for p in stack[:, y0 : y0 + window, x0 : x0 + window].astype(np.float64):
        p = np.nan_to_num(p)
        h = p - ndi.gaussian_filter(p, smooth)
        out.append(1.4826 * np.median(np.abs(h - np.median(h))) / p.std())
        sds.append(p.std())
    return float(np.median(out)), float(np.median(sds))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("clean")
    ap.add_argument("noisy")
    ap.add_argument("--like", nargs="+", required=True, help="measured projection stacks")
    ap.add_argument("--window", type=int, default=256)
    ap.add_argument("--smooth", type=float, default=1.5)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply the matched noise level")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--report-only", action="store_true")
    args = ap.parse_args()

    rd = [rho(tifffile.imread(f), args.window, args.smooth)[0] for f in args.like]
    for f, r in zip(args.like, rd):
        print(f"rho (data) {r:.4f}  {f}")
    clean = tifffile.imread(args.clean).astype(np.float64)
    r0, sd_win = rho(clean, args.window, args.smooth)
    rdata = float(np.mean(rd))
    # rho^2 (1 + s^2) = rho_0^2 + (W s)^2  (the std in the denominator contains the noise too)
    s = np.sqrt(max(rdata**2 - r0**2, 0.0) / (WHITE_RESPONSE**2 - rdata**2)) * args.scale
    print(f"rho (data, mean) {rdata:.4f}   rho (clean simulation) {r0:.4f}   "
          f"-> relative noise std s = {s:.4f} (x scale {args.scale:g})")
    if args.report_only:
        return
    rng = np.random.default_rng(args.seed)
    sigma = s * sd_win  # relative to the std over the central window, as in rho
    noisy = clean + rng.normal(0.0, sigma, clean.shape)
    tifffile.imwrite(args.noisy, noisy.astype(np.float32))
    print(f"sigma = {sigma:.4g} (window std {sd_win:.4g}); "
          f"rho (noisy) {rho(noisy, args.window, args.smooth)[0]:.4f} (target {rdata:.4f}); "
          f"wrote {args.noisy}")


if __name__ == "__main__":
    main()
