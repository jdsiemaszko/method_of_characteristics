from helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
import numpy as np

class GenericFlowElement():
    def __init__(self, v_plus, v_minus):
        self.v_plus = v_plus
        self.v_minus = v_minus

    @property
    def prandtl_meyer_angle(self):
        return (self.v_plus + self.v_minus) / 2

    @property
    def flow_direction(self):
        return (self.v_minus - self.v_plus) / 2

    @property
    def mach_number(self):
        return mach_from_prandtl_meyer(self.prandtl_meyer_angle)

    @property
    def mach_angle(self):
        return np.arcsin(1 / self.mach_number)

    @property
    def gamma_plus_direction(self):
        return self.flow_direction + self.mach_angle

    @property
    def gamma_minus_direction(self):
        return self.flow_direction - self.mach_angle

class FluidPoint(GenericFlowElement):
    def __init__(self, pos=None, v_plus=None, v_minus=None):
        self.pos = pos # x y coordinates

        # define point using two characteristics, if neither is provided, backpropagate!
        if v_plus is not None:
            self.v_plus = v_plus
        else:
            self.v_plus = self.backpropagate_gamma_plus()

        if v_minus is not None:
            self.v_minus = v_minus
        else:
            self.v_minus = self.backpropagate_gamma_minus()

        super().__init__(self.v_plus, self.v_minus)

    def backpropagate_gamma_plus(self):
        pass

    def backpropagate_gamma_minus(self):
        pass





