import numpy as np
from src.helper import prandtl_meyer_from_mach
from src.fluidPoint import GenericFlowElement, FluidPoint
from src.expansion_fan import JetExpansionFan
from src.cluster import GeometryCluster

import pathlib
import os

# ensure correct cwd
os.chdir(pathlib.Path(__file__).parent)
if "plots" not in os.listdir(os.getcwd()):
    os.makedirs(pathlib.Path(os.getcwd()).joinpath("plots"))

pm_angle_inlet = prandtl_meyer_from_mach(2.0, 1.4)
inlet_conditions = GenericFlowElement(pm_angle_inlet, pm_angle_inlet)

jef = JetExpansionFan(inlet=inlet_conditions, pressure_ratio=2.0, origin=(0, 1), NCHAR=4, gamma=1.4, type=-1)

inlet_points = [
    FluidPoint((0, yp), pm_angle_inlet, pm_angle_inlet, boundary="lower" if yp==0 else None) for yp in np.linspace(0, 1, 5, endpoint=False)
]


inlet_points.extend(jef.characteristic_origins)

gc = GeometryCluster(inlet_points)
gc.run(printFlag=True)
gc.plot_geometry(save=True)
