import numpy as np
from src.fluidPoint import FluidPoint

class Characteristic(): # a characteristic class
    def __init__(self, origin : FluidPoint, type):

        # a fluid point acting as the origin, defines the invariant, direction, and position of characteristic.
        self.origin = origin
        self.type = type  # 1 - gamma+, -1 - gamma-, 0 - flow direction

        self.gamma = origin.gamma
        self.ptot = origin.ptot

        match type:
            case 1:
                self.direction = origin.gamma_plus_direction
            case -1:
                self.direction = origin.gamma_minus_direction
            case 0:
                self.direction = origin.flow_direction

        self.end = None

        self.frontline_complement = None  # cache container

    @property
    def measure(self):  # flow direction distance of the characteristic from origin to end
        if self.end is None:
            return None
        else:
            return self.origin.flow_direction_dot_product(self.end)


    def update_bool_of_origin(self):
        match self.type:
            case 1:
                self.origin._gamma_plus_bool = True
            case -1:
                self.origin._gamma_minus_bool = True
            case 0:
                self.origin._gamma_zero_bool = True


    def __mul__(self, other) -> [FluidPoint, None]:
        # multiplication of two characteristics
        # defined as the intersection point!

        if self.origin * other.origin < 1e-10: # ignore coincident points!
            return None

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
        if np.abs(np.linalg.det(A) < 1e-10):
            return None  # No intersection, lines are parallel

        # Solve for t and u
        t, u = np.linalg.solve(A, b)

        # Compute the intersection point using the parametric equation of the first line
        position = tuple(np.array([x1, y1]) + t * dir1)

        if self.type == other.type: # intersection of characteristics of the same type
            return FluidPoint(position, None, None, gamma=self.gamma, ptot = self.ptot)

        # if self.type == 0 or other.type == 0: # if we're intersecting zero chars
        #     return FluidPoint(position, None, None)

        # 2) choose the correct invariants and initialize the intersection point
        if self.type == 1:
            vp = self.origin.v_plus
            if other.type == -1:
                vm = other.origin.v_minus
                b = None
            else: # zero char => reached upper boundary
                vm = -vp + 2 * other.origin.prandtl_meyer_angle
                b = "upper"
        elif self.type == -1:
            vm = self.origin.v_minus
            if other.type == 1:
                vp = other.origin.v_plus
                b = None
            else: # zero char => reached lower boundary
                vp = vm - 2 * other.origin.flow_direction
                b = "lower"
        else:
            if other.type == 1:  # zero char => reached upper boundary
                vp = other.origin.v_plus
                vm = -vp + 2 * self.origin.prandtl_meyer_angle
                b = "upper"
            elif other.type == -1: # zero char => reached lower boundary
                vm = other.origin.v_minus
                vp = vm - 2 * self.origin.flow_direction
                b = "lower"

        fp = FluidPoint(position, v_plus=vp, v_minus=vm, boundary=b, gamma=self.gamma, ptot = self.ptot)

        return fp


if __name__=="__main__":
    import matplotlib.pyplot as plt

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
    plt.show(save=True)