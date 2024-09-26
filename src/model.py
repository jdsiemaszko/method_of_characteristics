from fluidPoint import FluidPoint
from boundary import BoundaryLine
import numpy as np
from helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
import matplotlib.pyplot as plt
from characteristic import Characteristic

class Model():
    def __init__(self, config):



        self.N_CHARACTERISTICS = config["N_CHARACTERISTICS"]
        self.PA = config["pa"]
        self.PE = config["pe"]
        self.ME = config["Me"]
        self.PHI_E = config["phi_e"]
        self.GAMMA = config["gamma"]

        self.v_plus_e = prandtl_meyer_from_mach(self.ME, self.GAMMA) + self.PHI_E
        self.v_minus_e = self.v_plus_e - 2 * self.PHI_E

        self.entry_point = FluidPoint(None, self.v_plus_e, self.v_minus_e)

        self.point_A = (0, config["jet_width"])

        self.characteristics = [] # initialize characteristics as empty

    def run(self):

        self.expansion_fan(self.point_A)


    def expansion_fan(self, point):
        intermediate = self.PE / self.PA * (1 + self.GAMMA / 2 * self.ME**2)
        MA = np.sqrt((intermediate-1) * 2 / self.GAMMA)
        pm_A = prandtl_meyer_from_mach(MA, self.GAMMA)

        flow_direction_A = pm_A - self.v_plus_e
        region_A = FluidPoint(None, self.v_plus_e, pm_A + flow_direction_A)


        char_dir_min = - self.entry_point.mach_angle
        char_dir_max = -(region_A.mach_angle - region_A.flow_direction)

        # propagate characteristics

        for angle in np.linspace(char_dir_min, char_dir_max, self.N_CHARACTERISTICS):
            ch = Characteristic(self.point_A, angle, 0.0, 1)

            if self.characteristics and ch.direction > 0 and ch.direction > self.characteristics[-1].direction:
                break # stop shooting once the characteristics intersect somewhere upstream!
            self.characteristics.append(ch)

    def plot(self):
        fig, ax = plt.subplots()



