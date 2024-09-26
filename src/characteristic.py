
import numpy as np

class Characteristic():
    def __init__(self, origin, direction, invariant, type):
        self.origin = origin
        self.direction = direction
        self.invariant = invariant
        self.type = type # 0 - gamma+, 1 - gamma-

    def find_symmetry_point(self):
        x_sym = self.origin[0] - self.origin[1] / np.tan(self.direction)
        return (x_sym, 0.0)

    def find_reflection(self):
        new_origin = self.find_symmetry_point()
        new_direction = -self.direction
        new_invariant = 0.0
        new_type = not self.type

        return Characteristic(new_origin, new_direction, new_invariant, new_type)

    def __mul__(self, other):
        # multiplication of two characteristics
        # defined as the intersection point!

        delta_x = other.origin[0] - self.origin[0]
        delta_y = other.origin[1] - self.origin[1]

        tan1 = np.tan(self.direction)
        tan2 = np.tan(other.direction)


        dx1 = delta_x * tan2 / (tan1 - tan2)
        dy1 = tan1 * dx1

        return (self.origin[0] + dx1, self.origin[1] + dy1)