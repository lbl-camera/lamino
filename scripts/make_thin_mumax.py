#!/usr/bin/env python3
"""Thin version of a Mumax OVF: block-average along z and crop in-plane.

The simulation has 2 nm cells in z and 8 nm in x, y. The reconstruction code treats every
voxel as isotropic, so an emma-like sample (8-9 voxels of 8 nm) is made by averaging
--block consecutive z slices (default 4) into one voxel, --nblocks times, starting at slice
--z0. The in-plane crop [0:n] (default 319, odd, as m000000_odd.ovf) is kept. The averaged
m is what a coarser voxel holds, so |m| < 1 where m turns through the thickness (Neel caps).

Usage: make_thin_mumax.py src.ovf dst.ovf [--block 4] [--nblocks 9] [--z0 2] [--n 319]
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from wall_width import load_ovf  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--block", type=int, default=4)
    ap.add_argument("--nblocks", type=int, default=9)
    ap.add_argument("--z0", type=int, default=2)
    ap.add_argument("--n", type=int, default=319)
    args = ap.parse_args()

    m = load_ovf(args.src)  # (3, nz, ny, nx)
    nz0 = m.shape[1]
    z1 = args.z0 + args.block * args.nblocks
    if z1 > nz0:
        sys.exit(f"z range {args.z0}..{z1} exceeds the {nz0} slices of {args.src}")
    m = m[:, args.z0 : z1, : args.n, : args.n]
    m = m.reshape(3, args.nblocks, args.block, args.n, args.n).mean(axis=2)
    mag = np.linalg.norm(m, axis=0)
    print(f"thin volume {m.shape[1:]}, |m|: mean {mag.mean():.3f}, min {mag.min():.3f}")

    # header from the source with the new geometry
    with open(args.src) as f:
        head = f.read().partition("# Begin: Data Text")[0]
    dz = args.block * float(re.search(r"# zstepsize: (\S+)", head)[1])
    dx = float(re.search(r"# xstepsize: (\S+)", head)[1])
    new = {"xnodes": args.n, "ynodes": args.n, "znodes": args.nblocks,
           "xmax": args.n * dx, "ymax": args.n * dx, "zmax": args.nblocks * dz,
           "zstepsize": dz, "zbase": dz / 2}
    for k, v in new.items():
        head = re.sub(rf"# {k}: \S+", f"# {k}: {v:g}", head)
    flat = np.moveaxis(m, 0, -1).reshape(-1, 3)  # z, y, x order with x fastest
    with open(args.dst, "w") as f:
        f.write(head)
        f.write("# Begin: Data Text\n")
        np.savetxt(f, flat, fmt="%.8g %.8g %.8g ", newline="\n")
        f.write("# End: Data Text\n# End: Segment\n")
    print(f"wrote {args.dst}")


if __name__ == "__main__":
    main()
