import numpy as np
from src.helper import prandtl_meyer_from_mach
from src.fluidPoint import GenericFlowElement, FluidPoint
from src.expansionFan import JetExpansionFan
from src.cluster import GeometryCluster

import pathlib
import os

# ensure correct cwd
os.chdir(pathlib.Path(__file__).parent)
# create plots dir if not existing
if "plots" not in os.listdir(os.getcwd()):
    os.makedirs(pathlib.Path(os.getcwd()).joinpath("plots"))
# clean plots directory before running
plots_dir = os.path.join(os.curdir, 'plots')
[os.remove(os.path.join(plots_dir, f)) if os.path.isfile(os.path.join(plots_dir, f)) else os.rmdir(os.path.join(plots_dir, f)) for f in os.listdir(plots_dir)]

# INPUT CONDITONS
pressure_ratio = 2.0
jet_width = 1.0
gamma = 1.4
Mach_inlet = 2.5
atm_pressure = 101325 # Pa

# NUMERICS
N_fan = 20  # number of rays in the expansion fan
N_inlet = 20  # number of inlet points for characteristic propagation

pm_angle_inlet = prandtl_meyer_from_mach(Mach_inlet, gamma)
inlet_conditions = GenericFlowElement(pm_angle_inlet, pm_angle_inlet)

jef = JetExpansionFan(inlet=inlet_conditions, pressure_ratio=pressure_ratio, origin=(0, jet_width/2), NCHAR=N_fan, gamma=gamma, type=-1, pa=atm_pressure)
ptot = jef.total_pressure

inlet_points = [ # first point is a boundary!
    FluidPoint((0, yp), pm_angle_inlet, pm_angle_inlet, gamma=gamma, ptot = ptot,
               boundary="lower" if yp == 0 else None) for yp in np.linspace(0, jet_width/2, N_inlet, endpoint=False)
]

inlet_points.extend(jef.characteristic_origins)

gc = GeometryCluster(inlet_points)
gc.run(printFlag=True, plot_interval=20, max_iter=200, plotkwargs={
    'save' : True,
    'markers' : False,
    'plot_frontline' : True,
})

for attr in ['mach_number', 'pressure']:
    gc.plot_contours(attr, save=True, plot_characteristics=False, plot_frontline=True, plot_boundaries=True)

