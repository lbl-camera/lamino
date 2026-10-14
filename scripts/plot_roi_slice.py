"""
Plot an XZ-plane ROI slice from three per-component TIFF files.

The slice is normal to y (fixed y-index = cy).
The ROI in x is centred at cx with a configurable width (default: 4 × nz).
The full z-extent is always shown.

Component files are derived from the base path by appending 0/1/2.tiff:
    <base>0.tiff  – x-component
    <base>1.tiff  – y-component
    <base>2.tiff  – z-component

Examples:
    python plot_roi_slice.py
    python plot_roi_slice.py -f emma/665/recon -c 255 255 -x 105
    python plot_roi_slice.py -f emma/665/recon --arrows 20 --show
"""

import pathlib

import click
import matplotlib.pyplot as plt
import numpy as np
import tifffile


def _strip_tiff_ext(path: str) -> str:
    """Remove a trailing .tiff/.tif extension if present."""
    p = pathlib.Path(path)
    if p.suffix in {".tiff", ".tif"}:
        return str(p.with_suffix(""))
    return path


@click.command()
@click.option("-f", "--file", "base", default="emma/665/recon",
              show_default=True,
              help="Base path for component TIFFs (0/1/2.tiff appended).")
@click.option("-c", "--center", "center", nargs=2, type=int, default=(None, None),
              metavar="CX CY",
              help="ROI centre in voxel coordinates (default: volume centre).")
@click.option("-x", "--x-size", "x_size", type=int, default=None,
              help="ROI width in x (voxels). Default: 4 × nz.")
@click.option("--show", is_flag=True, default=False,
              help="Display the figure interactively in addition to saving.")
def main(base, center, x_size, show):
    # ----------------------------------------------------------------------- #
    # Load component TIFFs  →  each has shape (nz, ny, nx)
    # ----------------------------------------------------------------------- #
    base = _strip_tiff_ext(base)
    paths = [f"{base}{i}.tiff" for i in range(3)]

    comps = []
    for p in paths:
        click.echo(f"Reading {p}")
        comps.append(tifffile.imread(p).astype(np.float32))

    # Stack into (nz, ny, nx, 3)
    vecs = np.stack(comps, axis=-1)
    nz, ny, nx = vecs.shape[:3]

    # ----------------------------------------------------------------------- #
    # Resolve defaults that depend on volume size
    # ----------------------------------------------------------------------- #
    cx, cy = center
    if cx is None:
        cx = nx // 2
    if cy is None:
        cy = ny // 2
    if x_size is None:
        x_size = 4 * nz

    roi_half_x = x_size // 2
    x0, x1 = max(cx - roi_half_x, 0), min(cx + roi_half_x, nx)
    roi_width = x1 - x0

    click.echo(f"Volume dims   : {nx} x {ny} x {nz}")
    click.echo(f"ROI centre    : cx={cx},  y-slice={cy}")
    click.echo(f"ROI x-range   : [{x0}, {x1})  ({roi_width} voxels = {roi_width/nz:.1f}× z-thickness)")
    click.echo(f"Full z-range  : [0, {nz})")

    # ----------------------------------------------------------------------- #
    # Extract XZ ROI  →  shape (nz, roi_width, 3)
    # ----------------------------------------------------------------------- #
    roi = vecs[:, cy, x0:x1, :]                # [iz, ix_roi, component]
    mag = np.linalg.norm(roi, axis=-1)         # [iz, ix_roi]

    ix = np.arange(x0, x1)
    iz = np.arange(nz)
    extent = [ix[0], ix[-1], iz[0], iz[-1]]

    # ----------------------------------------------------------------------- #
    # Quiver subsampling
    # ----------------------------------------------------------------------- #
    step_x = 2
    step_z = 2

    qx = ix[::step_x]
    qz = iz[::step_z]
    qu = roi[::step_z, ::step_x, 0]            # x-component
    qw = roi[::step_z, ::step_x, 2]            # z-component
    QX, QZ = np.meshgrid(qx, qz)

    # ----------------------------------------------------------------------- #
    # Plot
    # ----------------------------------------------------------------------- #
    fig_w = max(6, roi_width / nz * 3)
    fig, ax = plt.subplots(figsize=(fig_w, 4))

    im = ax.imshow(
        mag,
        origin="lower",
        extent=extent,
        cmap="viridis",
        interpolation="bicubic",
        aspect="equal",
    )
    fig.colorbar(im, ax=ax, label="|magnetization|")

    ax.quiver(
        QX, QZ, qu, qw,
        color="black",
        pivot="mid",
        width=0.002,
        headwidth=3,
        headlength=4,
        alpha=0.85,
    )

    ax.set_title(
        f"XZ slice at y={cy}\n"
        f"ROI cx={cx}, x=[{x0},{x1}) ({roi_width} vox = {roi_width/nz:.1f}× nz)"
    )
    ax.set_xlabel("x (voxels)")
    ax.set_ylabel("z (voxels)")

    plt.tight_layout()
    out_png = f"roi_xz_cy{cy}_cx{cx}.png"
    plt.savefig(out_png, dpi=150)
    click.echo(f"Saved → {out_png}")

    if show:
        plt.show()


if __name__ == "__main__":
    main()
