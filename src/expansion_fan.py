
from src.fluidPoint import GenericFlowElement, FluidPoint
import numpy as np
from src.helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
from src.characteristic import Characteristic


class JetExpansionFan:
    def __init__(self, inlet : GenericFlowElement, pressure_ratio : float,
                 origin:tuple, NCHAR:int, gamma:float, type, pa = 1e6):

        self.inlet = inlet
        self.pressure_ratio = pressure_ratio
        self.pa = pa # atmospheric pressure (only passed for plotting purposes)
        self.gamma = gamma
        self.type = type
        self.origin = origin #tuple, not FluidPoint!
        self.compute_outlet() # computes fan outlet conditions
        self.nchar = NCHAR
        self.characteristics = [None] * NCHAR
        self.characteristic_origins = [None] * NCHAR
        self.initialize_characteristics()

        #  1 - downward facing (creating type 1 characteristics), -1 - upward facing (crating type -1 characteristics)
        self.type=type

    @property
    def total_pressure(self):
        return self.pa * (1 + self.gamma/2 * self.outlet.mach_number**2)

    def compute_outlet(self):
        intermediate = self.pressure_ratio * (1 + self.gamma / 2 * self.inlet.mach_number ** 2)
        mach_outlet = np.sqrt((intermediate - 1) * 2 / self.gamma)
        pm_outlet = prandtl_meyer_from_mach(mach_outlet, self.gamma)
        flow_direction_outlet = pm_outlet - self.inlet.v_plus

        self.outlet = GenericFlowElement(self.inlet.v_plus, pm_outlet + flow_direction_outlet)

    def initialize_characteristics(self):
        # assuming jet expansion increases flow angle
        char_dir_min = -(self.inlet.mach_angle - self.inlet.flow_direction)
        char_dir_max = -(self.outlet.mach_angle - self.outlet.flow_direction)

        for index, angle in enumerate(np.linspace(char_dir_min, char_dir_max, self.nchar)):

            # ASSUMING FLOW DIECTION CHANGES UNIFORMLY INSIDE THE FAN
            flow_direction_local = (angle - char_dir_min) / (char_dir_max - char_dir_min) *\
                                   (self.outlet.flow_direction - self.inlet.flow_direction) + self.inlet.flow_direction

            pm_local = flow_direction_local + self.inlet.v_plus
            v_minus_local = pm_local + flow_direction_local

            if index == self.nchar - 1: # last point is a boundary, all others are one-gamma only
                b = 'upper' if self.type==-1 else "lower"
            else:
                b = 'minus_only' if self.type==-1 else "plus_only"

            fp = FluidPoint(self.origin, v_plus=self.inlet.v_plus, v_minus=v_minus_local, boundary=b, gamma=self.gamma, ptot = self.total_pressure)
            char = Characteristic(fp, type=self.type) # gamma-
            self.characteristics[index] = char
            self.characteristic_origins[index] = fp
