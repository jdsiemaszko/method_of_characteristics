import matplotlib.pyplot as plt
import numpy as np
from characteristic import Characteristic
from fluidPoint import FluidPoint, GenericFlowElement

class GeometryCluster:
    def __init__(self, init_points):

        # init_points includes:
        # 1) arbitrary amount of inlet points
        # 2) arbitrary amount of expansion fan origin points (all coinciding, but with different invariants)

        # initialize point and char containers
        self.points = init_points
        self.gamma_plus = []
        self.gamma_minus = []

        self.jet_boundary = None # upper boundary
        self.symmetry_line = None # lower boundary

        self.frontline_points = init_points
        self.get_frontline_characteristics()

        self.closed_characteristics = []



    # def sort_frontline(self):
    #     # sort frontline, first by y position, then by flow direction
    #     sorted_frontline = sorted(self.frontline_points, key=lambda flp: (flp.pos[1], flp.flow_direction))
    #
    #     return list(sorted_frontline)

    def get_frontline_characteristics(self):
        self.front_gammma_plus = [None] * len(self.frontline_points)
        self.front_gammma_minus = [None] * len(self.frontline_points)

        for ind, p in enumerate(self.frontline_points):
            cp, cm = self.make_characteristics(p)
            self.front_gammma_plus[ind] = cp
            self.front_gammma_minus[ind] = cm

    def make_characteristics(self, point:FluidPoint):
        c_plus = Characteristic(point, type=0)
        c_minus = Characteristic(point, type=1)
        return c_plus, c_minus

    def find_first_intersection(self, char1, gamma_plus, gamma_minus):

        if char1.type: # if char1 is a gamma- type, search for intersects with gamma_plus
            gamma = gamma_plus
        else:
            gamma = gamma_minus
        intersects = [None] * len(gamma)
        for ind, g in enumerate(gamma):
            intersects[ind] = g * char1

        x0 = char1.origin.pos[0]

        # sort by x distance, reject any intersections left of char1.origin
        if all([inters.pos[0] <= x0 for inters in intersects]):
            # if all intersects are to the left, no new intersections found!
            return None

        sorted_intersect_and_gamma = sorted(zip(intersects, gamma), key=lambda ip, g: (ip.pos[0] if ip.pos[0] > x0 else 1e10))
        inter, g = sorted_intersect_and_gamma[0]
        return inter, g, char1

    def advance_frontline(self):
        # define new frontline points and store them in the cache

        new_frontline = []
        closest_intersects = [None] * len(self.frontline_points)
        for ind, p in enumerate(self.frontline_points):
            for char in self.make_characteristics(p):
                new_intersect, ch1, ch2 = self.find_first_intersection(char, self.front_gammma_plus, self.front_gammma_minus)
                if new_intersect is not None:
                    char.end = new_intersect
                    self.closed_characteristics.append(char)
                    new_frontline.append(new_intersect)


        # store new frontline and update front characteristics
        print('advanced frontline!')
        print('new frontline size: {}'.format(len(new_frontline)))
        self.frontline_points = new_frontline
        self.get_frontline_characteristics()
        self.points.extend(new_frontline)

    def plot_geometry(self, save = False):
        fig, ax = plt.subplots()

        ppos = [p.pos for p in self.points]

        ax.scatter(
            [pp[0] for pp in ppos],
            [pp[1] for pp in ppos],
            marker = 'x', color = 'k'
        )

        for closed_char in self.closed_characteristics:
            if closed_char.end is None:
                continue # if, for some reason, we reach a char that has no end point yet, skip it
            xplot = [closed_char.origin.pos[0], closed_char.end.pos[0]]
            yplot = [closed_char.origin.pos[1], closed_char.end.pos[1]]
            ax.plot(xplot, yplot, color='k', alpha = 0.3, linestyle='dashed')


        ax.set_ylim(0, 2)
        ax.set_xlim(0, max([p.pos[0] for p in self.frontline_points]))
        ax.grid()
        if save:
            plt.savefig(
                "/plots/{}.svg".format('geometry')
            )
        else:
            plt.show()



if __name__=="__main__":
    from expansion_fan import JetExpansionFan
    from helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
    pm_angle_inlet = prandtl_meyer_from_mach(2.0, 1.4)
    inlet_conditions = GenericFlowElement(pm_angle_inlet, pm_angle_inlet)

    jef = JetExpansionFan(inlet=inlet_conditions, pressure_ratio=2.0, origin=(0, 1), NCHAR=10, gamma=1.4, type=1)

    inlet_points = [
        FluidPoint((0, yp), pm_angle_inlet, pm_angle_inlet) for yp in np.linspace(0, 1, 10)
    ]

    inlet_points.extend(jef.characteristic_origins)

    gc = GeometryCluster(inlet_points)
    for i in range(3):
        gc.advance_frontline()
    gc.plot_geometry()


