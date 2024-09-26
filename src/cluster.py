from characteristic import Characteristic
from fluidPoint import FluidPoint

class GeometryCluster:
    def __init__(self, init_points, init_characteristics):

        self.points = init_points
        self.characteristics = init_characteristics
    def make_characteristics(self, point:FluidPoint):
        c_plus = Characteristic(point, type=0)
        c_minus = Characteristic(point, type=1)
        return c_plus, c_minus

    def next_characteristic_intersection(self, point):
        # find the next point to iterate
        pass


if __name__=="__main__":
    p1 = FluidPoint((0, 0), 0.5, 0.5)
    gc = GeometryCluster([p1], [])

    cp, cm = gc.make_characteristics(p1)
    p2 = cp * cm
    pass
