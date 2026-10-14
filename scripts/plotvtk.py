import sys
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

if __name__ == "__main__":
    # --- Read VTI file ---
    if len(sys.argv) > 1:
        vti_file = sys.argv[1]
    else:
        raise ValueError("Please provide a .vti file as a command-line argument.")

    out_file = sys.argv[2] if len(sys.argv) > 2 else "output.png"

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vti_file)
    reader.Update()

    # Rescale vectors so max magnitude = 1
    pd = reader.GetOutput().GetPointData()
    vecs = vtk_to_numpy(pd.GetArray("magetization"))  # shape (N, 3)
    magnitudes = np.linalg.norm(vecs, axis=1)
    max_mag = magnitudes.max()
    print(f"Max magnetization magnitude: {max_mag:.6f}")
    vecs_scaled = vecs / max_mag

    vtk_arr = numpy_to_vtk(vecs_scaled, deep=True)
    vtk_arr.SetName("magetization")
    pd.AddArray(vtk_arr)
    pd.SetActiveVectors("magetization")

    # --- Subsample points ---
    mask = vtk.vtkMaskPoints()
    mask.SetInputData(reader.GetOutput())
    mask.SetMaximumNumberOfPoints(3000)
    mask.RandomModeOn()
    mask.SetRandomModeType(1)  # uniform random sampling
    mask.Update()

    # --- Arrow glyph source ---
    arrow = vtk.vtkArrowSource()

    # --- Glyph filter ---
    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetInputConnection(mask.GetOutputPort())

    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    glyph.SetColorModeToColorByVector()
    glyph.OrientOn()
    glyph.SetScaleFactor(50)  # adjust as needed

    glyph.Update()

    glyph_out = glyph.GetOutput()
    print(f"Glyph output: {glyph_out.GetNumberOfPoints()} points, "
          f"{glyph_out.GetNumberOfCells()} cells")

    # --- Color by vector magnitude ---
    # Vectors are normalized so magnitude range is [0, 1]
    scalar_range = (0.0, 1.0)

    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.0)  # blue (low) -> red (high)
    lut.SetRange(scalar_range)
    lut.Build()

    # --- Mapper ---
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.SetLookupTable(lut)
    mapper.SetScalarRange(scalar_range)
    mapper.ScalarVisibilityOn()

    # --- Scalar bar (colormap legend) ---
    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetTitle("|magnetization|")
    scalarBar.SetNumberOfLabels(5)
    scalarBar.GetTitleTextProperty().SetColor(0, 0, 0)
    scalarBar.GetLabelTextProperty().SetColor(0, 0, 0)

    # --- Actor ---
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # --- Renderer ---
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddViewProp(scalarBar)
    renderer.SetBackground(1.0, 1.0, 1.0)

    # --- Render Window ---
    # Prefer EGL (Nvidia GPU, no display needed) > OSMesa (software) > default
    if hasattr(vtk, "vtkEGLRenderWindow"):
        renderWindow = vtk.vtkEGLRenderWindow()
    elif hasattr(vtk, "vtkOSOpenGLRenderWindow"):
        renderWindow = vtk.vtkOSOpenGLRenderWindow()
    else:
        renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(800, 600)

    # --- Render ---
    renderWindow.Render()

    # --- Capture image ---
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()

    # --- Save PNG ---
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(out_file)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()

    print(f"Saved visualization to {out_file}")
