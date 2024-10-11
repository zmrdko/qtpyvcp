import os
import numpy as np

from vtk import (
    vtkActor,
    vtkPoints,
    vtkPolyDataMapper,
    vtkPolyData,
    vtkDelaunay2D,
    vtkSmoothPolyDataFilter,
    vtkLinearSubdivisionFilter,
    vtkTransform
)

# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingContextOpenGL2
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkChartsCore import (
    vtkChartXYZ,
    vtkPlotSurface
)
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkFloatArray
from vtkmodules.vtkCommonDataModel import (
    vtkRectf,
    vtkTable,
    vtkVector2i
)
from vtkmodules.vtkRenderingContext2D import vtkContextMouseEvent
from vtkmodules.vtkViewsContext2D import vtkContextView

from vtkmodules.vtkFiltersCore import vtkSmoothPolyDataFilter
from vtkmodules.vtkFiltersModeling import vtkLoopSubdivisionFilter
from vtkmodules.vtkCommonCore import vtkLookupTable, vtkFloatArray
from vtkmodules.vtkRenderingCore import vtkPolyDataMapper


from qtpyvcp.utilities import logger

# from qtpyvcp.widgets.display_widgets.vtk_backplot.linuxcnc_datasource import LinuxCncDataSource
from qtpyvcp.utilities.settings import getSetting
from qtpyvcp.plugins import iterPlugins, getPlugin

LOG = logger.getLogger(__name__)


class PointsSurfaceActor(vtkActor):
    def __init__(self, datasource):
        super(PointsSurfaceActor, self).__init__()
        self.log = LOG
        self._datasource = datasource

        self.probe_results = []




        # self.axis = self._datasource.getAxis()
        # show_surface = getSetting('backplot.show-points-surface')
        # self.showSurface(show_surface and show_surface.value)

    def load_probe_results(self, filename):
        """Reads probe results from a text file and parses them into a list of X, Y, Z values.
        The first line should contain the current work coordinate offsets in the format: 'WCO X Y Z'."""
        if not os.path.exists(filename):
            self.log.error(f"File {filename} not found.")
            return  # Exit the method if the file doesn't exist

        try:
            with open(filename, 'r') as file:
                # Read the first line for offsets
                offsets_line = file.readline().strip()
                try:
                    # Split the line and extract values after 'WCO'
                    parts = offsets_line.split()
                    if parts[0] == "WCO" and len(parts) == 4:
                        # Convert the offset values to floats
                        self.mesh_x_offset = float(parts[1])
                        self.mesh_y_offset = float(parts[2])
                        self.mesh_z_offset = float(parts[3])
                        self.log.info(f"Loaded offsets: X={self.mesh_x_offset}, Y={self.mesh_y_offset}, Z={self.mesh_z_offset}")
                    else:
                        self.log.warning("Invalid format for offsets; expected 'WCO X Y Z'.")
                        return
                except ValueError as e:
                    self.log.error(f"Error parsing offsets: {e}")
                    return

                # Read the remaining lines for probe results
                line_count = 1  # Start counting lines from 1 since we've read the first line
                self.probe_results = []  # Reset the results

                for line in file:
                    line_count += 1
                    print(f"Raw line {line_count}: '{line.strip()}'")

                    row = list(map(float, line.strip().split()))
                    if len(row) == 3:
                        self.probe_results.append(row)
                    else:
                        self.log.warning(f"Invalid data format in line {line_count}: {line.strip()}")

                if not self.probe_results:
                    self.log.warning("No valid data found in the probe-results.txt file.")
                else:
                    self.log.info(f"Loaded {len(self.probe_results)} points from {filename}.")
                    print("Parsed probe results:")
                    for point in self.probe_results:
                        print(point)

        except FileNotFoundError:
            self.log.error(f"File {filename} not found.")
        except ValueError as e:
            self.log.error(f"Error parsing {filename}: {e}")

    def showSurface(self, show_surface):

        if show_surface:

            self.VisibilityOn()
            self.log.info("SHOW POINTS SURFACE")
            self.load_probe_results('probe-results.txt')

            # Check if probe_results is not empty
            if not self.probe_results:
                self.log.error("No probe results available to show.")
                return

            # Transform the probe results to a structured 2D array
            structured_results = np.array(self.probe_results)

            num_points_x = structured_results.shape[0]

            # Extract X, Y, Z coordinates
            x_coords = structured_results[:, 0]
            y_coords = structured_results[:, 1]
            z_coords = structured_results[:, 2]

            x_min, x_max = x_coords.min(), x_coords.max()
            y_min, y_max = y_coords.min(), y_coords.max()

            # Create a vtkPoints object and add the points
            vtk_points = vtkPoints()
            for point in structured_results:
                vtk_points.InsertNextPoint(point)

            # Create a vtkPolyData object and set the points
            polydata = vtkPolyData()
            polydata.SetPoints(vtk_points)

            # **Assign Z-values as scalar data** for coloring
            z_scalars = vtkFloatArray()
            z_scalars.SetName("Z-Scalars")
            for z in z_coords:
                z_scalars.InsertNextValue(z)

            polydata.GetPointData().SetScalars(z_scalars)  # Set Z as scalar data

            # Create a vtkDelaunay2D filter to generate the mesh
            delaunay = vtkDelaunay2D()
            delaunay.SetInputData(polydata)
            delaunay.Update()

            # Optional: Apply a subdivision filter for higher-resolution mesh
            subdivision_filter = vtkLinearSubdivisionFilter()
            subdivision_filter.SetInputConnection(delaunay.GetOutputPort())
            subdivision_filter.SetNumberOfSubdivisions(3)  # Adjust the number of subdivisions for smoother surface
            subdivision_filter.Update()

            # Apply smoothing filter (bicubic-like)
            smooth_filter = vtkSmoothPolyDataFilter()
            smooth_filter.SetInputConnection(subdivision_filter.GetOutputPort())  # Connect to subdivision filter
            smooth_filter.SetNumberOfIterations(100)  # Increase for smoother results
            smooth_filter.SetRelaxationFactor(0.5)  # Lower factor for smoother but less aggressive smoothing
            smooth_filter.Update()

            # Get min and max Z for color mapping (heatmap)
            z_min, z_max = z_coords.min(), z_coords.max()

            # Create a color lookup table (LUT) for heatmap colors
            lut = vtkLookupTable()
            lut.SetTableRange(z_min, z_max)
            lut.SetHueRange(0.0, 0.67)  # From red to blue (heatmap style)
            lut.Build()

            # Create a mapper and set the input connection to the smoothed surface
            mapper = vtkPolyDataMapper()
            mapper.SetInputConnection(smooth_filter.GetOutputPort())
            mapper.SetLookupTable(lut)
            mapper.SetScalarRange(z_min, z_max)
            mapper.ScalarVisibilityOn()

            # Apply mesh transformations
            x_offset = self.mesh_x_offset
            y_offset = self.mesh_y_offset
            z_offset = self.mesh_z_offset

            mesh_transform = vtkTransform()
            mesh_transform.Translate(x_offset, y_offset, z_offset)

            # Apply the transformations and set the mapper
            self.SetUserTransform(mesh_transform)
            self.SetMapper(mapper)
            self.GetProperty().SetPointSize(10)
            sm_opacity = 1-getSetting('backplot.surface-map-transparency').value/100
            print(f"opacity: {sm_opacity}")
            self.GetProperty().SetOpacity(sm_opacity)  # Set transparency (0.0 = fully transparent, 1.0 = fully opaque)

        else:
            self.VisibilityOff()
            self.log.info("HIDE POINTS SURFACE")



    # def run(self):
    #         self.view.GetInteractor().Start()