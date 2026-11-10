import numpy as np
import matplotlib.pyplot as plt
import tifffile
from pathlib import Path


def plot_vectors(vectors, output_path, islice=None, roi=None, axis=None):
    """
    Plots the vectors as arrows and saves the plot to the specified output path.
    """

    if len(vectors.shape) != 3:
        raise ValueError("input is not a tupe of 3D arrays")

    # Determine axis to plot
    if axis is None:
        axis = 0  # Default to the first axis

    if islice is None:
        islice = vectors[0].shape[axis] // 2 

    # slice the vectors along the specified axis
    if axis == 0:
        mx = vectors[0][islice, :, :]
        my = vectors[1][islice, :, :]
        mz = vectors[2][islice, :, :]
    elif axis == 1:
        mx = vectors[0][:, islice, :]
        my = vectors[1][:, islice, :]
        mz = vectors[2][:, islice, :]
    elif axis == 2:
        mx = vectors[0][:, :, islice]
        my = vectors[1][:, :, islice]
        mz = vectors[2][:, :, islice]
    else:
        raise ValueError("Invalid axis. Must be 0, 1, or 2.")

    # magnitude of the vectors
    magnitude = np.sqrt(mx**2 + my**2 + mz**2)


def load_vectors_from_tiff(path: pathlib.Path):
    """
    Loads the vectors from a TIFF file and returns them as a tuple of 3D arrays.
    """
    vectors = tifffile.imread(str(path))
    return vectors
