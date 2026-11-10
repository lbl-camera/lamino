"""
Check consistency between TIFF and VTI reconstructions.

TIFFs store individual vector components in (z, y, x) axis order,
named as <basename>0.tiff, <basename>1.tiff, <basename>2.tiff.
VTI file is named <basename>.vti.

Usage:
    python check_consistency.py <path/to/basename>

Example:
    python check_consistency.py emma/665/recon
"""

import sys
import numpy as np
import tifffile
import vtk
from vtk.util.numpy_support import vtk_to_numpy

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <path/to/basename>")
    sys.exit(1)

basename = sys.argv[1]

# --------------------------------------------------------------------------- #
# Load TIFFs — each is (nz, ny, nx)
# --------------------------------------------------------------------------- #
tiff_components = []
for i in range(3):
    path = f"{basename}{i}.tiff"
    arr = tifffile.imread(path)
    tiff_components.append(arr)
    print(f"  {path}: shape={arr.shape}, dtype={arr.dtype}")

tiff_stack = np.stack(tiff_components, axis=-1)   # (nz, ny, nx, 3)
nz, ny, nx, _ = tiff_stack.shape
print(f"\nTIFF stack (z,y,x,comp): {tiff_stack.shape}")

# --------------------------------------------------------------------------- #
# Load VTI
# --------------------------------------------------------------------------- #
vti_path = f"{basename}.vti"
reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(vti_path)
reader.Update()
image = reader.GetOutput()

dims = image.GetDimensions()   # (nx, ny, nz) in VTK convention
print(f"VTI dimensions (VTK x,y,z): {dims}")

pd = image.GetPointData()
vtk_arr = pd.GetArray(0)
print(f"VTI array: name='{vtk_arr.GetName()}', "
      f"components={vtk_arr.GetNumberOfComponents()}, "
      f"tuples={vtk_arr.GetNumberOfTuples()}")

# vtk_to_numpy gives a flat (ntuples, ncomp) array; VTK iterates x-fastest
# Reshape to (nx, ny, nz, 3) in Fortran order, then transpose to (nz, ny, nx, 3)
vti_flat = vtk_to_numpy(vtk_arr).reshape(dims[0], dims[1], dims[2], 3, order='F')
vti_zyx  = np.transpose(vti_flat, (2, 1, 0, 3))   # (nz, ny, nx, 3)
print(f"VTI array reshaped (z,y,x,comp): {vti_zyx.shape}\n")

# --------------------------------------------------------------------------- #
# Consistency checks
# --------------------------------------------------------------------------- #
assert tiff_stack.shape == vti_zyx.shape, (
    f"Shape mismatch: TIFF {tiff_stack.shape} vs VTI {vti_zyx.shape}")

component_names = ["x", "y", "z"]
all_ok = True

for c, cname in enumerate(component_names):
    t = tiff_stack[..., c]
    v = vti_zyx[..., c]

    abs_diff = np.abs(t - v)
    max_diff  = abs_diff.max()
    mean_diff = abs_diff.mean()
    rel_max   = max_diff / (np.abs(t).max() + 1e-30)

    match = np.allclose(t, v, rtol=1e-5, atol=1e-7)
    status = "PASS" if match else "FAIL"
    all_ok = all_ok and match

    print(f"Component {cname} [{status}]")
    print(f"  max |diff|  = {max_diff:.3e}   (rel {rel_max:.3e})")
    print(f"  mean |diff| = {mean_diff:.3e}")
    print(f"  TIFF range  [{t.min():.6f}, {t.max():.6f}]")
    print(f"  VTI  range  [{v.min():.6f}, {v.max():.6f}]")

print()
if all_ok:
    print("All components match — TIFF and VTI are consistent.")
else:
    print("Mismatch detected — TIFF and VTI are NOT fully consistent.")
