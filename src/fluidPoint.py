from src.helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
import numpy as np

class GenericFlowElement():
    def __init__(self, v_plus, v_minus, gamma=1.4, ptot = 1e6):
        self.v_plus = v_plus
        self.v_minus = v_minus
        self.gamma = gamma
        self.ptot = ptot

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

    @property
    def pressure_over_total_pressure(self):
        return 1 / (1 + self.gamma / 2 * self.mach_number**2)

    @property
    def pressure(self):
        return self.pressure_over_total_pressure * self.ptot

class FluidPoint(GenericFlowElement): # generic flow element with position added!
    def __init__(self, pos, v_plus=None, v_minus=None, boundary=None, gamma=1.4, ptot = 1e6):

        self.pos = pos # x y coordinates

        self.boundary = boundary # flag to check if point is on the boundary

        # booleans checking if characteristics have been tried
        self._gamma_plus_bool = False
        self._gamma_minus_bool = False
        self._gamma_zero_bool = False

        # match the characteristics to shoot based on the boundary condition
        if self.boundary == "upper":
            self._gamma_plus_bool = True
        if self.boundary == "lower":
            self._gamma_minus_bool = True
        if self.boundary == "plus_only":
            self._gamma_zero_bool = True
            self._gamma_minus_bool = True
        if self.boundary == "minus_only":
            self._gamma_zero_bool = True
            self._gamma_plus_bool = True
        if self.boundary is None:
            self._gamma_zero_bool = True

        # accepted values: "lower", "upper"

        self.ending_characteristics = set({}) # store reference to all characteristics ending at this point

        super().__init__(v_plus, v_minus, gamma, ptot)

    @property
    def all_chars_exhausted(self):
        return all((self._gamma_plus_bool, self._gamma_minus_bool, self._gamma_zero_bool))

    def flow_direction_dot_product(self, other):
        flow_direction_vec = np.array([
            np.cos(self.flow_direction),
            np.sin(self.flow_direction)
        ])
        delta_vec = np.array([
            other.pos[0] - self.pos[0],
            other.pos[1] - self.pos[1],
        ])

        return np.dot(flow_direction_vec, delta_vec)

    def __mul__(self, other) -> float:
        # multiplication of points => return distance squared

        return (self.pos[0]-other.pos[0])**2 + (self.pos[1] - other.pos[1])**2


