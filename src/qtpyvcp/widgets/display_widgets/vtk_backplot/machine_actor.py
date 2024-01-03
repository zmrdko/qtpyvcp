
import math

from pprint import pprint
from _collections import defaultdict

import vtk.qt
from vtk.util.colors import *

from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersSources import vtkLineSource
from vtkmodules.vtkFiltersSources import vtkCylinderSource
from vtkmodules.vtkRenderingCore import vtkPolyDataMapper

from qtpyvcp.utilities import logger
from collections import OrderedDict

LOG = logger.getLogger(__name__)


class MachineActor(vtk.vtkCubeAxesActor):

    def __init__(self, linuxcncDataSource):
        super(MachineActor, self).__init__()
        self._datasource = linuxcncDataSource
        axis = self._datasource.getAxis()
        units = self._datasource.getProgramUnits()

        x_max = axis[0]["max_position_limit"]
        x_min = axis[0]["min_position_limit"]

        y_max = axis[1]["max_position_limit"]
        y_min = axis[1]["min_position_limit"]

        z_max = axis[2]["max_position_limit"]
        z_min = axis[2]["min_position_limit"]

        self.SetBounds(x_min, x_max, y_min, y_max, z_min, z_max)

        self.SetXLabelFormat("%6.3f")
        self.SetYLabelFormat("%6.3f")
        self.SetZLabelFormat("%6.3f")

        self.SetFlyModeToStaticEdges()

        self.GetTitleTextProperty(0).SetColor(1.0, 0.0, 0.0)
        self.GetLabelTextProperty(0).SetColor(1.0, 0.0, 0.0)

        self.GetTitleTextProperty(1).SetColor(0.0, 1.0, 0.0)
        self.GetLabelTextProperty(1).SetColor(0.0, 1.0, 0.0)

        self.GetTitleTextProperty(2).SetColor(0.0, 0.0, 1.0)
        self.GetLabelTextProperty(2).SetColor(0.0, 0.0, 1.0)

        self.SetXUnits(units)
        self.SetYUnits(units)
        self.SetZUnits(units)

        self.DrawXGridlinesOn()
        self.DrawYGridlinesOn()
        self.DrawZGridlinesOn()

        self.SetGridLineLocation(self.VTK_GRID_LINES_FURTHEST)

        self.GetXAxesGridlinesProperty().SetColor(0.0, 0.0, 0.0)
        self.GetYAxesGridlinesProperty().SetColor(0.0, 0.0, 0.0)
        self.GetZAxesGridlinesProperty().SetColor(0.0, 0.0, 0.0)

    def showMachineTicks(self, ticks):
        if ticks:
            self.XAxisTickVisibilityOn()
            self.YAxisTickVisibilityOn()
            self.ZAxisTickVisibilityOn()
        else:
            self.XAxisTickVisibilityOff()
            self.YAxisTickVisibilityOff()
            self.ZAxisTickVisibilityOff()

    def showMachineBounds(self, bounds):
        if bounds:
            self.XAxisVisibilityOn()
            self.YAxisVisibilityOn()
            self.ZAxisVisibilityOn()
        else:
            self.XAxisVisibilityOff()
            self.YAxisVisibilityOff()
            self.ZAxisVisibilityOff()

    def showMachineLabels(self, labels):
        if labels:
            self.XAxisLabelVisibilityOn()
            self.YAxisLabelVisibilityOn()
            self.ZAxisLabelVisibilityOn()
        else:
            self.XAxisLabelVisibilityOff()
            self.YAxisLabelVisibilityOff()
            self.ZAxisLabelVisibilityOff()

    def showGridlines(self, grid):
        if grid:
            self.DrawXGridlinesOn()
            self.DrawYGridlinesOn()
            self.DrawZGridlinesOn()
        else:
            self.DrawXGridlinesOff()
            self.DrawYGridlinesOff()
            self.DrawZGridlinesOff()


class MachinePart(vtk.vtkAssembly):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.part_axis = None
        self.part_type = None

    def SetPartAxis(self, attr):
        self.part_axis = attr
        
    def SetPartType(self, attr):
        self.part_type = attr

    def GetPartAxis(self):
        return self.part_axis
    
    def GetPartType(self):
        return self.part_type


class MachinePartsASM(vtk.vtkAssembly):
    
    def __init__(self, parts):
        super(MachinePartsASM, self).__init__()
        
        print("########")
        print("NEW Machine")
        
        # print(f"{parts}")
        previous_asm = None
        
        parts_dict = OrderedDict()
        previous_depth = 0
        branch_num = 0
        
        # for depth, part_root, part_data in self.items_recursive(parts):
        for part in self.items_recursive(parts, self):
            print(f"{part}")
            
        print(f"vtkAssembly: {self}")


    def create_part(self, data):
        
        part_id = data.get("id")
        part_model = data.get("model")
        part_type = data.get("type")
        part_position = data.get("position")
        part_origin = data.get("origin")
        part_axis = data.get("axis")
        part_joint = data.get("joint")
        
        print(f"part_id:\t\t{part_id}")
        print(f"part_model:\t\t{part_model}")
        print(f"part_type:\t\t{part_type}")
        print(f"part_position:\t\t{part_position}")
        print(f"part_origin:\t\t{part_origin}")
        print(f"part_axis:\t\t{part_axis}")
        print(f"part_joint:\t\t{part_joint}")
        
        part_source = vtk.vtkSTLReader()
        part_source.SetFileName(part_model)
        part_source.Update()
        
        part_mapper = vtk.vtkPolyDataMapper()
        part_mapper.SetInputConnection(part_source.GetOutputPort())
        
        part_actor = vtk.vtkActor()
        
        part_actor.SetMapper(part_mapper)
        
        part_actor.GetProperty().SetColor(1, 0, 1)
        part_actor.GetProperty().SetDiffuseColor(0.9, 0.9, 0.9)
        part_actor.GetProperty().SetDiffuse(.8)
        part_actor.GetProperty().SetSpecular(.5)
        part_actor.GetProperty().SetSpecularColor(1.0, 1.0, 1.0)
        part_actor.GetProperty().SetSpecularPower(30.0)
        
        part_actor.SetPosition(part_position[0], part_position[1], part_position[2])
        part_actor.SetOrigin(part_origin[0], part_origin[1], part_origin[2])
        
        tmp_assembly = MachinePart()
        tmp_assembly.SetPartAxis(part_axis)
        tmp_assembly.SetPartType(part_type)
        
        tmp_assembly.AddPart(part_actor)
        
        return tmp_assembly

    def items_recursive(self, d, parent):
        
        for _, v in d.items():
            if isinstance(v, dict):
                tmp_part = self.create_part(v)
                parent.AddPart(tmp_part)
                yield tmp_part
                for p in self.items_recursive(v, tmp_part):
                    parent.AddPart(tmp_part)
                    yield tmp_part
        
                