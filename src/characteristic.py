
import numpy as np
from fluidPoint import GenericFlowElement, FluidPoint

class Characteristic(): # a characteristic class
    def __init__(self, origin : FluidPoint, type):

        # a fluid point acting as the origin, defines the invariant, direction, and position of characteristic.
        self.origin = origin
        self.type = type  # 0 - gamma+, 1 - gamma-
        self.direction = origin.gamma_plus_direction if not type else origin.gamma_minus_direction
        self.end = None

    def find_symmetry_point(self, y_reflect=0):
        x_sym = self.origin.pos[0] - (self.origin.pos[1]-y_reflect) / np.tan(self.direction)

        flow_direction_sym = 0 # symmetry condition!
        pm_angle_sym = self.origin.v_plus + flow_direction_sym if not self.type\
            else self.origin.v_minus - flow_direction_sym

        vp= pm_angle_sym - flow_direction_sym
        vm = pm_angle_sym + flow_direction_sym

        pt = FluidPoint((x_sym, 0), vp, vm)
        return pt

    def reflection_over_symmetry(self, y_reflect = 0):
        new_origin = self.find_symmetry_point(y_reflect)
        new_type = not self.type

        return Characteristic(new_origin, new_type)

    def __mul__(self, other) -> [FluidPoint, None]:
        # multiplication of two characteristics
        # defined as the intersection point!

        # 0) if characteristics have the same direction and different origin => no intersection
        if self.direction == other.direction:
            return None

        # 1) Compute the position of the intersection point

        delta_x = other.origin.pos[0] - self.origin.pos[0]
        delta_y = other.origin.pos[1] - self.origin.pos[1]

        tan1 = np.tan(self.direction)
        tan2 = np.tan(other.direction)

        dx1 = (delta_x * tan2 - delta_y) / (tan1 - tan2)
        dy1 = tan1 * dx1

        position = (self.origin.pos[0] + dx1, self.origin.pos[1] + dy1)

        if self.type == other.type: # intersection of characteristics of the same type
            return FluidPoint(position, None, None)

        # 2) choose the correct invariants and initialize the intersection point
        vp = self.origin.v_plus if not self.type else other.origin.v_plus
        vm = self.origin.v_minus if self.type else other.origin.v_minus

        fp = FluidPoint(position, v_plus=vp, v_minus=vm)

        # 3) store the end point (convenient for plotting)
        self.end = fp
        other.end = fp

        return fp


if __name__=="__main__":
    fp1 = FluidPoint((0, 0), 0.5, 0.3)
    c1 = Characteristic(fp1, type=0)

    fp2 = FluidPoint((0, 1), 0.7, 0.6)
    c2 = Characteristic(fp2, type=1)

    fp3 = c1 * c2
