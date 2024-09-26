
import numpy as np
from fluidPoint import GenericFlowElement, FluidPoint

class Characteristic(): # a characteristic class
    def __init__(self, origin : FluidPoint, type):

        # a fluid point acting as the origin, defines the invariant, direction, and position of characteristic.
        self.origin = origin
        self.type = type  # 0 - gamma+, 1 - gamma-
        self.direction = origin.gamma_plus_direction if not type else origin.gamma_minus_direction

    # def find_symmetry_point(self):
    #     x_sym = self.origin.pos[0] - self.origin.pos[1] / np.tan(self.direction)
    #
    #     vp = self.origin.v_plus # reflection only works for gamma +
    #     vm = 0 # ???
    #     pt = FluidPoint((x_sym, 0), vp, vm)
    #     return pt
    #
    # def find_reflection(self):
    #     new_origin = self.find_symmetry_point()
    #     new_direction = -self.direction
    #     new_type = not self.type
    #
    #     return Characteristic(new_origin, new_type)

    def __mul__(self, other) -> [FluidPoint, None]:
        # multiplication of two characteristics
        # defined as the intersection point!

        # 0) if characteristics are of the same type, return None (model no longer valid)
        if self.type == other.type:
            return None

        # 1) Compute the position of the intersection point

        delta_x = other.origin.pos[0] - self.origin.pos[0]
        delta_y = other.origin.pos[1] - self.origin.pos[1]

        tan1 = np.tan(self.direction)
        tan2 = np.tan(other.direction)

        dx1 = (delta_x * tan2 - delta_y) / (tan1 - tan2)
        dy1 = tan1 * dx1

        position = (self.origin.pos[0] + dx1, self.origin.pos[1] + dy1)

        # 2) pass the correct invariants to the intersection point
        vp = self.origin.v_plus if not self.type else other.origin.v_plus
        vm = self.origin.v_minus if self.type else other.origin.v_minus
        return FluidPoint(position, v_plus=vp, v_minus=vm)


if __name__=="__main__":
    fp1 = FluidPoint((0, 0), 0.5, 0.3)
    c1 = Characteristic(fp1, type=0)

    fp2 = FluidPoint((0, 1), 0.7, 0.6)
    c2 = Characteristic(fp2, type=1)

    fp3 = c1 * c2

    pass
