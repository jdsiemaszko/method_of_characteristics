
import numpy as np
from fluidPoint import GenericFlowElement, FluidPoint

class Characteristic():
    def __init__(self, origin : FluidPoint, direction : float, type):
        self.origin = origin
        self.direction = direction
        self.type = type  # 0 - gamma+, 1 - gamma-

    def find_symmetry_point(self):
        x_sym = self.origin.pos[0] - self.origin.pos[1] / np.tan(self.direction)

        vp = self.origin.v_plus # reflection only works for gamma +
        vm = 0 # ???
        pt = FluidPoint((x_sym, 0), vp, vm)
        return pt

    def find_reflection(self):
        new_origin = self.find_symmetry_point()
        new_direction = -self.direction
        new_type = not self.type

        return Characteristic(new_origin, new_direction, new_type)

    def __mul__(self, other):
        # multiplication of two characteristics
        # defined as the intersection point!

        delta_x = other.origin[0] - self.origin.pos[0]
        delta_y = other.origin[1] - self.origin.pos[1]

        tan1 = np.tan(self.direction)
        tan2 = np.tan(other.direction)


        dx1 = (delta_x * tan2 - delta_y) / (tan1 - tan2)
        dy1 = tan1 * dx1

        position = (self.origin.pos[0] + dx1, self.origin.pos[1] + dy1)

        if self.type!=other.type: # intersection of + and - characteristics => valid intersection point
            vp = self.origin.v_plus if not self.type else other.origin.v_plus
            vm = self.origin.v_minus if self.type else other.origin.v_minus
            return FluidPoint(position, v_plus=vp, v_minus=vm)

        else: # intersection of two + or two - characteristics => shock formation! model invalid!
            return None

