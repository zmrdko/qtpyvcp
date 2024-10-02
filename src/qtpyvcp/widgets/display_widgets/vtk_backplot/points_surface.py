import os
import numpy as np

from vtk import (
    vtkActor,
    vtkPoints,
    vtkPolyDataMapper,
    vtkPolyData,
    vtkDelaunay2D,
    vtkSmoothPolyDataFilter,
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
        self.mesh_x_offset = getSetting("backplot.mesh-x-offset").value  # Initialize mesh_x_offset
        self.mesh_y_offset = getSetting("backplot.mesh-y-offset").value  # Initialize mesh_y_offset
        self.mesh_z_offset = getSetting("backplot.mesh-z-offset").value  # Initialize mesh_z_offset

        self.axis = self._datasource.getAxis()
        # show_surface = getSetting('backplot.show-points-surface')
        # self.showSurface(show_surface and show_surface.value)

    def load_probe_results(self, filename):
        """Reads probe results from a text file and parses them into a list of X, Y, Z values."""
        try:
            with open(filename, 'r') as file:
                line_count = 0
                total_lines = sum(1 for line in file)  # Count total lines
                file.seek(0)  # Reset file pointer to the beginning

                print(f"Total lines in file: {total_lines}")

                for line in file:
                    line_count += 1
                    # Print the raw line to debug
                    print(f"Raw line {line_count}: '{line.strip()}'")

                    # Strip any leading/trailing whitespace and split by spaces
                    row = list(map(float, line.strip().split()))
                    
                    # Check the parsed row length
                    if len(row) == 3:
                        self.probe_results.append(row)
                    else:
                        self.log.warning(f"Invalid data format in line {line_count}: {line.strip()}")

            if not self.probe_results:
                self.log.warning("No valid data found in the probe-results.txt file.")
            else:
                self.log.info(f"Loaded {len(self.probe_results)} points from {filename}.")
                # Output the parsed results to verify them
                print("Parsed probe results:")
                for point in self.probe_results:
                    print(point)  # Each point is a list of [X, Y, Z]

        except FileNotFoundError:
            self.log.error(f"File {filename} not found.")
        except ValueError as e:
            self.log.error(f"Error parsing {filename}: {e}")

    def showSurface(self, show_surface):

        if show_surface:
            self.log.info("SHOW POINTS SURFACE ")
            self.load_probe_results('probe-results.txt')

            # Check if probe_results is not empty
            if not self.probe_results:
                self.log.error("No probe results available to show.")
                return
            
            # Ensure probe_results is structured correctly
            # Transform the probe results to a structured 2D array
            structured_results = np.array(self.probe_results)
            
            num_points_x = structured_results.shape[0]
            num_points_y = structured_results.shape[1]

            # Log the sizes
            self.log.info(f"Number of points in X: {num_points_x}, Y: {num_points_y}")

            # Extract X, Y, Z coordinates
            x_coords = structured_results[:, 0]
            y_coords = structured_results[:, 1]
            z_coords = structured_results[:, 2]

            x_min, x_max = x_coords.min(), x_coords.max()
            y_min, y_max = y_coords.min(), y_coords.max()

            # Create a grid of X and Y based on the actual coordinates
            x_range = np.linspace(x_min, x_max, num_points_x)
            y_range = np.linspace(y_min, y_max, num_points_y)

            # Initialize the points array
            points_array = np.zeros((num_points_x * num_points_y, 3))

            # Generate the points with elevation and translated coordinates
            for idx, (x, y, z) in enumerate(structured_results):
                points_array[idx] = [x, y, z]

            # Print the points array for verification
            print("Points Array:")
            print(points_array)

            # Create a vtkPoints object and add the points
            vtk_points = vtkPoints()
            for point in points_array:
                vtk_points.InsertNextPoint(point)

            # Create a vtkPolyData object and set the points
            polydata = vtkPolyData()
            polydata.SetPoints(vtk_points)

            # Create a vtkDelaunay2D filter to generate the mesh
            delaunay = vtkDelaunay2D()
            delaunay.SetInputData(polydata)
            delaunay.Update()

            # Apply mesh transformations
            x_offset = self.mesh_x_offset+(x_max-x_min)/2
            y_offset = self.mesh_y_offset+(y_max-y_min)/2
            z_offset = self.mesh_z_offset

            mesh_transform = vtkTransform()
            mesh_transform.Translate(x_offset, y_offset, z_offset)

            # Create a mapper and set the input connection
            mapper = vtkPolyDataMapper()
            mapper.SetInputConnection(delaunay.GetOutputPort())

            self.SetUserTransform(mesh_transform)
            self.SetMapper(mapper)
            self.GetProperty().SetPointSize(10)

        else:
            self.log.info("HIDE POINTS SURFACE ")



    # def run(self):
    #         self.view.GetInteractor().Start()