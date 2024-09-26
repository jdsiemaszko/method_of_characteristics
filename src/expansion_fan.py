
from fluidPoint import GenericFlowElement, FluidPoint
import numpy as np
from helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
from characteristic import Characteristic


class JetExpansionFan:
    def __int__(self, inlet : GenericFlowElement, pressure_ratio : float, origin:tuple, NCHAR:int, gamma:float, type):

        self.inlet = inlet
        self.pressure_ratio = pressure_ratio
        self.gamma = gamma
        self.origin = origin
        self.compute_outlet() # computes fan outlet conditions
        self.nchar = NCHAR
        self.characteristics = np.empty(NCHAR, dtype = Characteristic)
        self.initialize_characteristics()

        #  0 - downward facing (creating type 0 characteristics), 1 - upward facing (crating type 1 characteristics)
        self.type=type
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

            fp = FluidPoint(self.origin, v_plus=self.inlet.v_plus, v_minus=v_minus_local)
            char =  Characteristic(fp, angle, type=self.type) # gamma-
            self.characteristics[index] = char