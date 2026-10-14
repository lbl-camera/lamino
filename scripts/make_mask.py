"""
Make a binary support mask from the central slice of a scalar reconstruction.

The mask is 1 inside a centred --region square and 0 outside it and inside the
fiducial holes. The fiducials are circles found in the central slice with a
Hough transform (canny edges). The mask is written as a float32 0/1 TIFF of
size --shape (the volume size recon_dims[1:]) that tomocam-lamino reads via
`mask = "..."` in [recon_params]. A PNG preview is written next to it.

Examples:
    python make_mask.py ../laminocam/emma/665/tv/recon.tiff emma/665/mask.tiff
    python make_mask.py recon.tiff mask.tiff --region 399 399 --n-holes 4
"""


import pathlib

import click
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.feature import canny
from skimage.transform import hough_circle, hough_circle_peaks


def centre_crop(img: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Crop img to shape about its centre (same placement as crop3d)."""
    h, w = shape
    if h > img.shape[0] or w > img.shape[1]:
        raise click.ClickException(
            f"slice {img.shape} is smaller than --shape {shape}")
    r0 = (img.shape[0] - h) // 2
    c0 = (img.shape[1] - w) // 2
    return img[r0:r0 + h, c0:c0 + w]


def find_fiducials(img, n_holes, rmin, rmax):
    """Return (cy, cx, r) of the n_holes strongest circles in img."""
    scale = np.abs(img).max()
    edges = canny(img / scale, sigma=2, low_threshold=0.02, high_threshold=0.08)
    radii = np.arange(rmin, rmax + 1)
    acc = hough_circle(edges, radii)
    _, cx, cy, r = hough_circle_peaks(acc, radii, total_num_peaks=n_holes,
                                      min_xdistance=2 * rmin,
                                      min_ydistance=2 * rmin)
    return list(zip(cy, cx, r))


@click.command()
@click.argument("recon", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(dir_okay=False))
@click.option("--shape", nargs=2, type=int, default=(511, 511),
              show_default=True, metavar="H W",
              help="Output size, must equal recon_dims[1:] of the config.")
@click.option("--region", nargs=2, type=int, default=(399, 399),
              show_default=True, metavar="H W",
              help="Centred region inside which the mask is 1.")
@click.option("--slice", "islice", type=int, default=None,
              help="Slice used to find the fiducials (default: central).")
@click.option("--n-holes", type=int, default=4, show_default=True,
              help="Number of fiducial holes.")
@click.option("--radius", nargs=2, type=int, default=(18, 32),
              show_default=True, metavar="MIN MAX",
              help="Fiducial radius search range in pixels.")
@click.option("--margin", type=float, default=3.0, show_default=True,
              help="Extra radius added to each hole (covers the rim).")
def main(recon, output, shape, region, islice, n_holes, radius, margin):
    vol = np.asarray(tifffile.imread(recon), dtype=np.float32)
    if vol.ndim == 2:
        vol = vol[None]
    k = vol.shape[0] // 2 if islice is None else islice
    img = centre_crop(vol[k], tuple(shape))

    mask = np.zeros(shape, dtype=bool)
    r0 = (shape[0] - region[0]) // 2
    c0 = (shape[1] - region[1]) // 2
    mask[r0:r0 + region[0], c0:c0 + region[1]] = True

    yy, xx = np.mgrid[:shape[0], :shape[1]]
    holes = find_fiducials(img, n_holes, *radius)
    for cy, cx, r in holes:
        mask &= (yy - cy) ** 2 + (xx - cx) ** 2 > (r + margin) ** 2
        print(f"fiducial at (row {cy}, col {cx}), radius {r}")

    out = pathlib.Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out, mask.astype(np.float32))

    fig, ax = plt.subplots(1, 2, figsize=(9, 4.5))
    ax[0].imshow(img, cmap="gray")
    ax[0].set_title(f"recon, slice {k}")
    ax[1].imshow(mask, cmap="gray")
    ax[1].set_title("mask")
    for a in ax:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=150)
    print(f"mask fraction {mask.mean():.3f}, shape {mask.shape} -> {out}")


if __name__ == "__main__":
    main()
