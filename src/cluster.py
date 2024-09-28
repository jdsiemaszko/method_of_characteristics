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

        self.jet_boundary = None # upper boundary
        self.symmetry_line = None # lower boundary

        self.frontline_points = init_points # init points should be an iterable (ideally a set)
        self.get_frontline_characteristics()

        self.dead_points = set({})
        self.dead_characteristics = set({})



    # def sort_frontline(self):
    #     # sort frontline, first by y position, then by flow direction
    #     sorted_frontline = sorted(self.frontline_points, key=lambda flp: (flp.pos[1], flp.flow_direction))
    #
    #     return list(sorted_frontline)

    def get_frontline_characteristics(self):

        self.front_gammma_plus = []
        self.front_gammma_minus = []
        self.front_boundary = []

        for ind, p in enumerate(self.frontline_points):
            cp, cm, c0 = self.make_characteristics(p)
            if cp is not None:
                self.front_gammma_plus.append(cp)
            if cm is not None:
                self.front_gammma_minus.append(cm)
            if c0 is not None:
                self.front_boundary.append(c0)

        self.frontline_characteristics = self.front_gammma_plus + self.front_gammma_minus + self.front_boundary

    def make_characteristics(self, point : FluidPoint):

        # if point.boundary == "upper":
        #     c_plus = None
        #     c_minus = Characteristic(point, type=-1)
        #     c_0 = Characteristic(point, type=0)
        # elif point.boundary == "lower":
        #     c_plus = Characteristic(point, type=1)
        #     c_minus = None
        #     c_0 = Characteristic(point, type=0)
        # else:
        #     c_plus = Characteristic(point, type=1)
        #     c_minus = Characteristic(point, type=-1)
        #     c_0 = None

        if point._gamma_plus_bool:
            c_plus = None
        else:
            c_plus = Characteristic(point, type=1)
        if point._gamma_minus_bool:
            c_minus = None
        else:
            c_minus = Characteristic(point, type=-1)
        if point._gamma_zero_bool:
            c_0 = None
        else:
            c_0 = Characteristic(point, type=0)
        # c_plus = Characteristic(point, type=1)
        # c_minus = Characteristic(point, type=-1)
        # c_0 = Characteristic(point, type=0)
        return c_plus, c_minus, c_0

    def find_first_intersection(self, char1):

        # if char1.type==1: # if char1 is a gamma+ type, search for intersects with gamma- and gamma0
        #     gamma = gamma_plus + gamma_0
        # elif char1.type == -1:
        #     gamma = gamma_minus + gamma_0
        # else:
        #
        gamma = self.frontline_characteristics
        gamma_valid = []
        intersects = []
        for ind, g in enumerate(gamma):
            inter = g * char1
            if inter is not None and inter.v_plus is not None and\
                    char1.origin.flow_direction_dot_product(inter) > 0:
                    # and\
                    #     inter.flow_direction_dot_product(g.origin) > 0: # if the intersection is valid and downstream for both characteristics
                intersects.append(g * char1)
                gamma_valid.append(g)

        x0 = char1.origin.pos[0]



        #if no such intersects exist, return nones
        if not intersects: # empty list is false-y
            return None, None, None

        # sort by flow_distance!

        sorted_intersect_and_gamma = sorted(zip(intersects, gamma_valid), key=lambda tuple: (char1.origin.flow_direction_dot_product(tuple[0])))
        inter, g = sorted_intersect_and_gamma[0]

        return inter, g, char1

    def advance_frontline(self):
        # define new frontline points and store them in the cache

        for char in self.frontline_characteristics:
            new_intersect, ch_other, _ = self.find_first_intersection(char)

            # NOTE: new intersect an ch_other could be None!
            if new_intersect is not None and new_intersect.v_plus is not None:
                # check if this inmtersect is indeed the closest for both characteristics
                if ch_other.measure is None or ch_other.origin.flow_direction_dot_product(new_intersect) < ch_other.measure:
                    ch_other.end = new_intersect
                    ch_other.frontline_complement = char
                    new_intersect.ending_characteristics.append(ch_other)

                # this should ALWAYS be true (check)
                # if char.length_squared is None or new_intersect * char.origin < char.length_squared:

                char.end = new_intersect
                char.frontline_complement = ch_other
                new_intersect.ending_characteristics.append(char)



        new_frontline_points = set({}) # a set of points!
        new_dead_points = set({}) # a set of points!
        new_dead_characteristics = set({}) # a set of chars!

        for char in self.frontline_characteristics:
            # if __frontline_complement relation is mutual (end point is closest to both points)
            if char.frontline_complement is not None and char.frontline_complement.frontline_complement == char:
                new_frontline_points.add(char.end) # add end point to new front
                char.update_bool_of_origin()
                new_dead_characteristics.add(char)
            # elif char.frontline_complement is not None and char.frontline_complement.frontline_complement != char: # frontline complement is not mutual
            #     new_frontline_points.add(char.origin) # add origin point to new front, in hopes of finding a mutual complement in the next advance
            #
            # else:  # if no complement found at all don't add anything to the new front!
            #     new_dead_points.add(char.origin) # origin is now dead as well

        for point in self.frontline_points: # add any non-dead points in the current frontline to the new one
            if not point.all_chars_exhausted:
                new_frontline_points.add(point)
            else:
                new_dead_points.add(point) # kill the point
                for dead_char in point.ending_characteristics:
                    new_dead_characteristics.add(dead_char)


        # store new frontline and update front characteristics
        print('advanced frontline!')
        print('new frontline size: {}'.format(len(new_frontline_points)))
        self.frontline_points = new_frontline_points
        self.get_frontline_characteristics()
        self.dead_points = self.dead_points.union(new_dead_points)
        self.dead_characteristics = self.dead_characteristics.union(new_dead_characteristics)

    def plot_geometry(self, save = False):
        fig, ax = plt.subplots()

        # for p in list(self.dead_points) + list(self.frontline_points):
        #
        #     ax.text(p.pos[0], p.pos[1], str(p.boundary))


        for dead_char in self.dead_characteristics:
            if dead_char.end is None:
                continue # if, for some reason, we reach a char that has no end point yet, skip it

            xplot = [dead_char.origin.pos[0], dead_char.end.pos[0]]
            yplot = [dead_char.origin.pos[1], dead_char.end.pos[1]]
            ax.plot(xplot, yplot, color='k', alpha = 0.5, linestyle='dashed', marker='x')

        for fl_char in self.frontline_characteristics:
            xplot = [fl_char.origin.pos[0]]
            yplot = [fl_char.origin.pos[1]]
            ax.plot(xplot, yplot, color='r', alpha=0.5, linestyle='dashed', marker='x')


        # ax.set_ylim(0, 2)
        ax.set_xlim(0, max([p.pos[0] for p in self.frontline_points]))

        ax.set_xlabel('x')
        ax.set_xlabel('y')

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

    jef = JetExpansionFan(inlet=inlet_conditions, pressure_ratio=2.0, origin=(0, 1), NCHAR=10, gamma=1.4, type=-1)

    inlet_points = [
        FluidPoint((0, yp), pm_angle_inlet, pm_angle_inlet, boundary="lower" if yp==0 else None) for yp in np.linspace(0, 1, 5, endpoint=False)
    ]


    inlet_points.extend(jef.characteristic_origins)

    gc = GeometryCluster(inlet_points)
    for i in range(50):
        gc.advance_frontline()
    gc.plot_geometry()


