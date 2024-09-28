import matplotlib.pyplot as plt
import numpy as np
from fluidPoint import GenericFlowElement, FluidPoint

class Characteristic(): # a characteristic class
    def __init__(self, origin : FluidPoint, type):

        # a fluid point acting as the origin, defines the invariant, direction, and position of characteristic.
        self.origin = origin
        self.type = type  # 1 - gamma+, -1 - gamma-, 0 - flow direction

        match type:
            case 1:
                self.direction = origin.gamma_plus_direction
            case -1:
                self.direction = origin.gamma_minus_direction
            case 0:
                self.direction = origin.flow_direction

        self.end = None

        self.frontline_complement = None  # cache container

    def find_symmetry_point(self, y_reflect=0):
        if self.type == 0:
            return None

        x_sym = self.origin.pos[0] - (self.origin.pos[1]-y_reflect) / np.tan(self.direction)

        flow_direction_sym = 0 # symmetry condition!
        if self.type==1:
            pm_angle_sym = self.origin.v_plus + flow_direction_sym
        elif self.type==-1:
            pm_angle_sym = self.origin.v_minus - flow_direction_sym

        vp = pm_angle_sym - flow_direction_sym
        vm = pm_angle_sym + flow_direction_sym

        pt = FluidPoint((x_sym, 0), vp, vm)
        return pt

    def reflection_over_symmetry(self, y_reflect = 0):
        new_origin = self.find_symmetry_point(y_reflect)
        new_type = -self.type

        return Characteristic(new_origin, new_type)

    @property
    def length_squared(self):
        if self.end is None:
            return None
        else:
            return self.end * self.origin # distance squared between end and origin

    def __mul__(self, other) -> [FluidPoint, None]:
        # multiplication of two characteristics
        # defined as the intersection point!

        # Unpack origin points
        x1, y1 = self.origin.pos
        x2, y2 = other.origin.pos

        theta1 = self.direction
        theta2 = other.direction

        # Direction vectors for both lines
        dir1 = np.array([np.cos(theta1), np.sin(theta1)])
        dir2 = np.array([np.cos(theta2), np.sin(theta2)])

        # Set up the system of linear equations A * [t, u]^T = b
        A = np.array([[dir1[0], -dir2[0]],
                      [dir1[1], -dir2[1]]])

        b = np.array([x2 - x1, y2 - y1])

        # Check if the determinant is zero (lines are parallel)
        if np.linalg.det(A) == 0:
            return None  # No intersection, lines are parallel

        # Solve for t and u
        t, u = np.linalg.solve(A, b)

        # Compute the intersection point using the parametric equation of the first line
        position = tuple(np.array([x1, y1]) + t * dir1)

        if self.type == other.type: # intersection of characteristics of the same type
            return FluidPoint(position, None, None)

        # if self.type == 0 or other.type == 0: # if we're intersecting zero chars
        #     return FluidPoint(position, None, None)

        # 2) choose the correct invariants and initialize the intersection point
        if self.type == 1:
            vp = self.origin.v_plus
            if other.type == -1:
                vm = other.origin.v_minus
            else: # zero char
                vm = vp + 2 * other.origin.flow_direction
        elif self.type == -1:
            vm = self.origin.v_minus
            if other.type == 1:
                vp = other.origin.v_plus
            else: # zero char
                vp = vm - 2 * other.origin.flow_direction
        else: # self is a zero char
            if other.type == 1:
                vp = other.origin.v_plus
                vm = vp + 2 * self.origin.flow_direction
            elif other.type == -1:
                vm = other.origin.v_minus
                vp = vm - 2 * self.origin.flow_direction

        fp = FluidPoint(position, v_plus=vp, v_minus=vm)

        return fp


if __name__=="__main__":
    fp1 = FluidPoint((0, 0), 0.5, 0.3)
    c1 = Characteristic(fp1, type=1)

    fp2 = FluidPoint((0.5, 1), 0.7, 0.6)
    c2 = Characteristic(fp2, type=-1)

    fp3 = c1 * c2

    fig, ax = plt.subplots()

    ax.set_aspect('equal')
    ax.grid()
    print(np.rad2deg(fp1.gamma_plus_direction))
    print(np.rad2deg(fp2.gamma_minus_direction))

    ax.scatter(*fp1.pos)
    ax.scatter(*fp2.pos)
    ax.scatter(*fp3.pos)
    plt.show()