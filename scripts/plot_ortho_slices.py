

import pathlib
import click
import matplotlib.pyplot as plt
import numpy as np
import tifffile


def _strip_tiff_ext(path: str) -> str:
    p = pathlib.Path(path)
    if p.suffix in {".tiff", ".tif"}:
        return str(p.with_suffix(""))
    return path


@click.command()
@click.argument("tiffs", type=str, default="recon")
    
@click.option("--axis", type=str, default='z', show_default=True,
              help="Normal axis for slice")
@click.option("--slice", "islc", type=int, default=None, show_default=True,
                help="Slice index for XY plane (default: 0 for 'z', ny//2 for 'y', nx//2 for 'x').")
@click.option("-o", "--output", "out_png", default=None,
              help="Output PNG path (default: ortho_islc<axis>.png).")
@click.option("--show", is_flag=True, default=False,
              help="Display the figure interactively in addition to saving.")
def main(tiffs, axis, islc, out_png, show):
    # ----------------------------------------------------------------------- #
    # Load TIFFs  →  stack to (nz, ny, nx, 3)
    # ----------------------------------------------------------------------- #
    base = _strip_tiff_ext(tiffs)
    paths = [f"{base}{i}.tiff" for i in range(3)]

    components = []
    for p in paths:
        click.echo(f"Reading {p}")
        components.append(tifffile.imread(p).astype(np.float32))

    vecs = np.stack(components, axis=-1)   # (nz, ny, nx, 3)
    nz, ny, nx = vecs.shape[:3]

    # ----------------------------------------------------------------------- #
    # Resolve slice indices
    # ----------------------------------------------------------------------- #
    if axis == 'z':
        if islc is None:
            islc = 0
        elif islc >= nz:
            click.echo(f"Warning: iz={islc} exceeds max index {nz-1}, clamping to {nz-1}.")
            islc = nz - 1
    elif axis == 'y':
        if islc is None:
            islc = ny // 2
        elif islc >= ny:
            click.echo(f"Warning: iy={islc} exceeds max index {ny-1}, clamping to {ny-1}.")
            islc = ny - 1
    else:
        if islc is None:
            islc = nx // 2
        elif islc >= nx:
            click.echo(f"Warning: ix={islc} exceeds max index {nx-1}, clamping to {nx-1}.")
            islc = nx - 1

    # clamp
    click.echo(f"Volume dims : {nx} x {ny} x {nz}")
    if axis == 'z':
        click.echo(f"XY plane    : z={islc}")
    elif axis == 'y':
        click.echo(f"XZ plane    : y={islc}")
    else:
        click.echo(f"YZ plane    : x={islc}")

    ix = np.arange(nx)
    iy = np.arange(ny)
    iz = np.arange(nz)

    # ----------------------------------------------------------------------- #
    # Extract slices and compute magnitude
    # ----------------------------------------------------------------------- #
    if axis == 'z':
        s = vecs[islc, :, :, :].squeeze()  # (ny, nx, 3)
    elif axis == 'y':
        s = vecs[:, islc, :, :].squeeze()  # (nz, nx, 3)
    else:
        s = vecs[:, :, islc, :].squeeze()  # (nz, ny, 3)

    # normalize for consistent color scaling across panels
    mag = np.linalg.norm(s, axis=-1)  # (nr, nc)
    mag = mag / np.max(mag)

    # ----------------------------------------------------------------------- #
    # Plot – XY on left (full height), XZ/YZ stacked on right
    # ----------------------------------------------------------------------- #
    # adust figure size based on slice dimensions
    fig, axs = plt.subplots(1, 1)
    axs.set_aspect('equal')
    axs.imshow(mag, origin='lower', cmap='viridis')
    axs.set_xlabel("X (voxels)")
    axs.set_ylabel("Y (voxels)")
    plt.tight_layout()

    if out_png is None:
        out_png = f"ortho_{axis}{islc}.png"
    plt.savefig(out_png, dpi=300)
    click.echo(f"Saved → {out_png}")

    if show:
        plt.show()


if __name__ == "__main__":
    main()
