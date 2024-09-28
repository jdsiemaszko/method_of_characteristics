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

class FluidPoint(GenericFlowElement): # generic flow element with position added!
    def __init__(self, pos, v_plus=None, v_minus=None, boundary=None):

        self.pos = pos # x y coordinates

        self.v_plus = v_plus
        self.v_minus = v_minus

        self.boundary = boundary # flag to check if point is on the boundary

        # booleans checking if characteristics have been tried
        self._gamma_plus_bool = False
        self._gamma_minus_bool = False
        # accepted values: "lower", "upper"

        super().__init__(self.v_plus, self.v_minus)

    def __mul__(self, other) -> float:
        # multiplication of points => return distance squared

        return (self.pos[0]-other.pos[0])**2 + (self.pos[1] - other.pos[1])**2


