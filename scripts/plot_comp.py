import sys
import tifffile
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python plot_comp.py <path>")
        sys.exit(1)

    data_path = sys.argv[1]
    if not Path(data_path).is_dir():
        print(f"Error: {data_path} is not a directory.")
        sys.exit(1)
    xcomp_path = Path(data_path) / "recon0.tiff"
    ycomp_path = Path(data_path) / "recon1.tiff"
    zcomp_path = Path(data_path) / "recon2.tiff"

    xcomp = tifffile.imread(xcomp_path)
    ycomp = tifffile.imread(ycomp_path)
    zcomp = tifffile.imread(zcomp_path)

    # normalize the magnitude of the components
    mag = np.sqrt(xcomp**2 + ycomp**2 + zcomp**2)
    mag[mag == 0] = 1  # avoid division by zero
    max_mag = np.max(mag)
    if max_mag > 0:
        xcomp /= max_mag
        ycomp /= max_mag
        zcomp /= max_mag

    # set add  dictionary named imgs to basedir of the input file paths
    pngdir = Path(xcomp_path).parent / "png"
    Path(pngdir).mkdir(parents=True, exist_ok=True)

    xcomp_png = pngdir / "xcomp.png"
    ycomp_png = pngdir / "ycomp.png"
    zcomp_png = pngdir / "zcomp.png"

    # sum over the slow axis, and collapse the 3 components into a single image
    xcomp2d = np.sum(xcomp, axis=0)
    ycomp2d = np.sum(ycomp, axis=0)
    zcomp2d = np.sum(zcomp, axis=0)

    # save individual components as png files
    plt.imsave(xcomp_png, xcomp2d, cmap='gray')
    plt.imsave(ycomp_png, ycomp2d, cmap='gray')
    plt.imsave(zcomp_png, zcomp2d, cmap='gray')

    """"
    # plot the 3 components side by side
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    axs[0].imshow(xcomp2d, cmap='gray')
    axs[0].set_title('X Component')
    axs[0].axis('off')
    axs[1].imshow(ycomp2d, cmap='gray')
    axs[1].set_title('Y Component')
    axs[1].axis('off')
    axs[2].imshow(zcomp2d, cmap='gray')
    axs[2].set_title('Z Component')
    axs[2].axis('off')
    plt.tight_layout()
    plt.show()
    """
