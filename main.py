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



# INPUT CONDITONS
pressure_ratio = 2.0
jet_width = 1.0
gamma = 1.4
Mach_inlet = 2.0

# NUMERICS
N_fan = 10
N_inlet = 5

pm_angle_inlet = prandtl_meyer_from_mach(Mach_inlet, gamma)
inlet_conditions = GenericFlowElement(pm_angle_inlet, pm_angle_inlet)

jef = JetExpansionFan(inlet=inlet_conditions, pressure_ratio=pressure_ratio, origin=(0, jet_width), NCHAR=N_fan, gamma=gamma, type=-1)

inlet_points = [ # first point is a boundary!
    FluidPoint((0, yp), pm_angle_inlet, pm_angle_inlet,
               boundary="lower" if yp == 0 else None) for yp in np.linspace(0, jet_width, N_inlet, endpoint=False)
]

inlet_points.extend(jef.characteristic_origins)

gc = GeometryCluster(inlet_points)
gc.run(printFlag=True, plot_interval=10, max_iter=200, plotkwargs={
    'save' : True,
    'markers' : False
})
gc.plot_geometry(save=True, markers=False)
