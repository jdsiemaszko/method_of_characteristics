import os

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from src.characteristic import Characteristic
from src.fluidPoint import FluidPoint, GenericFlowElement

class GeometryCluster:
    def __init__(self, init_points):

        self.frontline_points = init_points # init points should be an iterable (ideally a set)
        self.frontline_characteristics = []
        self.get_frontline_characteristics()

        self.dead_points = set({})
        self.dead_characteristics = set({})

        self.iter = 0

    def get_frontline_characteristics(self):

        front_gammma_plus = []
        front_gammma_minus = []
        front_boundary = []

        for ind, p in enumerate(self.frontline_points):
            cp, cm, c0 = self.make_characteristics(p)
            if cp is not None:
                front_gammma_plus.append(cp)
            if cm is not None:
                front_gammma_minus.append(cm)
            if c0 is not None:
                front_boundary.append(c0)

        self.frontline_characteristics = front_gammma_plus + front_gammma_minus + front_boundary
        # sort frontline by x value!
        self.frontline_characteristics = sorted(self.frontline_characteristics, key = lambda x: x.origin.pos[0])

    def make_characteristics(self, point : FluidPoint):
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
        return c_plus, c_minus, c_0

    def find_first_intersection(self, char1):

        gamma = self.frontline_characteristics
        gamma_valid = []
        intersects = []
        for ind, g in enumerate(gamma):
            inter = g * char1
            if inter is not None and\
                    char1.origin.flow_direction_dot_product(inter) > 0\
                and\
                    g.origin.flow_direction_dot_product(inter) > 0:
            # make sure:
            # 1) intersection exists
            # 1) we don't backtrack from char1.origin
            # 3) we don't backtrack from g.origin.origin (if any)
            #     if not g.origin.ending_characteristics or g.origin.ending_characteristics.pop().origin.flow_direction_dot_product(inter) > 0:

                intersects.append(g * char1)
                gamma_valid.append(g)

        #if no such intersects exist, return nones
        if not intersects: # empty list is false-y
            return None, None, None

        # sort by flow_distance!

        sorted_intersect_and_gamma = sorted(zip(intersects, gamma_valid),
                                            key=lambda tuple: (char1.origin.flow_direction_dot_product(tuple[0])))
        inter, g = sorted_intersect_and_gamma[0]

        return inter, g, char1

    def find_first_dead_intersection(self, char1):

        gamma = self.dead_characteristics
        gamma_valid = []
        intersects = []
        for ind, g in enumerate(gamma):
            inter= g * char1
            if g.end.pos[0] > char1.origin.pos[0]: # don't need to bother with g to the left of char1.origin
                if inter is not None and g.origin != char1.origin and\
                        char1.origin.flow_direction_dot_product(inter) > 0 and\
                        g.origin.flow_direction_dot_product(inter) > 0 and\
                        g.end.flow_direction_dot_product(inter) < 0:

                # make sure:
                # 1) intersection exists
                # 1) we don't backtrack from char1.origin
                # 3) we don't backtrack from g.origin
                # 4) we DO backtrack from g.end (g is dead after all)
                # 5) g.origin is not same as char1.origin (trivial case)

                    intersects.append(g * char1)
                    gamma_valid.append(g)

        #if no such intersects exist, return nones
        if not intersects: # empty list is false-y
            return None, None, None

        # sort by flow_distance!

        sorted_intersect_and_gamma = sorted(zip(intersects, gamma_valid),
                                            key=lambda tuple: (char1.origin.flow_direction_dot_product(tuple[0])))
        inter, g = sorted_intersect_and_gamma[0]

        return inter, g, char1

    def advance_frontline(self, printFlag = False):
        # define new frontline points and store them in the cache

        new_frontline_points = set({}) # a set of points!
        new_dead_points = set({}) # a set of points!
        new_dead_characteristics = set({}) # a set of chars!

        stopFlag = False

        for char in self.frontline_characteristics:
            new_intersect, ch_other, _ = self.find_first_intersection(char)

            # NOTE: new intersect and ch_other could be None!
            if new_intersect is not None:

                # check if this inmtersect is indeed the closest for both characteristics
                if ch_other.measure is None or\
                        ch_other.origin.flow_direction_dot_product(new_intersect) < ch_other.measure:

                    ch_other.end = new_intersect
                    ch_other.frontline_complement = char

                # this should ALWAYS be true (check)
                # if char.length_squared is None or new_intersect * char.origin < char.length_squared:

                char.end = new_intersect
                char.frontline_complement = ch_other

        for char in self.frontline_characteristics:
            # if __frontline_complement relation is mutual (end point is closest to both points)
            if char.frontline_complement is not None and char.frontline_complement.frontline_complement == char:
                if char.end.v_plus is not None: # if there is no shock at char.end, continue the frontline
                    new_frontline_points.add(char.end) # add end point to new front
                else: # if we detect a shock, add end to dead points (cannot continue the model)
                    new_dead_points.add(char.end)

                    if printFlag:
                        print('shockwave formation detected at ({:.2f}, {:.2f})'.format(char.end.pos[0],
                                                                                        char.end.pos[1]))
                    stopFlag = True

                new_dead_characteristics.add(char)
                char.update_bool_of_origin()
                char.end.ending_characteristics.add(char)

        for point in self.frontline_points:
            if not point.all_chars_exhausted: # if we didn't find new intersections in the above loop
                if printFlag:
                    print('failed to find frontline intersect for point ({:.2f}, {:.2f}), checking dead characteristics!'.format(point.pos[0], point.pos[1]))
                for char in self.make_characteristics(point):
                    if char is not None:  # only check valid chars
                        inter, char2, _ = self.find_first_dead_intersection(char) # last chance - check the dead chars
                        if inter is None: # we didn't find any intersections with dead chars! => try again in the next frontline
                            new_frontline_points.add(point)
                        else: # we did find an intersection with dead chars!
                            if inter.v_plus is None: # but it's a shock!
                                new_dead_points.add(inter)
                                if printFlag:
                                    print('shockwave formation detected at ({:.2f}, {:.2f})'.format(inter.pos[0],
                                                                                                    inter.pos[1]))
                                stopFlag = True
                            else: # inter is a valid intersection point => add it to the frontline and make sure we don't shoot char2 again!
                                match char2.type:
                                    case 1:
                                        inter._gamma_plus_bool = True
                                    case -1:
                                        inter._gamma_minus_bool = True
                                    case 0:
                                        inter._gamma_zero_bool = True

                                new_frontline_points.add(inter)


            else:  # kill the point
                new_dead_points.add(point)
                for dead_char in point.ending_characteristics:
                    new_dead_characteristics.add(dead_char)

        # store new frontline and update front characteristics
        if printFlag:
            print('new frontline size: {}'.format(len(new_frontline_points)))

        self.frontline_points = new_frontline_points
        if stopFlag: # if the frontline is empty
            if printFlag:
                print('stopping model')
            return True

        self.get_frontline_characteristics()
        self.dead_points = self.dead_points.union(new_dead_points)
        self.dead_characteristics = self.dead_characteristics.union(new_dead_characteristics)

        return False  # continue the run function

    def run(self, max_iter = 100, printFlag = False, plot_interval=0, plotkwargs={'save' : True, 'markers' : False}):
        breakLoop = False
        while self.iter < max_iter and not breakLoop:
            print('iteration {}: advancing frontline points'.format(self.iter))
            print('current reach: {:.2f}'.format(
                max([p.pos[0] for p in self.frontline_points])
            ))
            breakLoop = self.advance_frontline(printFlag=printFlag)
            self.iter +=1
            if plot_interval > 0 and self.iter % plot_interval == 0:
                if printFlag:
                    print('plotting current geometry')
                self.plot_geometry(**plotkwargs)
        if plot_interval > 0:
            if printFlag:
                print('plotting current geometry')
            self.plot_geometry(**plotkwargs)


    def plot_geometry(self, save = False, markers=True, plot_frontline=True):
        fig, ax = plt.subplots(figsize = (8, 6))

        for dead_char in self.dead_characteristics:

            if dead_char.end is None:
                continue # if, for some reason, we reach a char that has no end point yet, skip it

            xplot = [dead_char.origin.pos[0], dead_char.end.pos[0]]
            yplot = [dead_char.origin.pos[1], dead_char.end.pos[1]]

            if dead_char.type in [1, -1]: # normal chars colored in grey, other configs can be passed through the advance function
                ax.plot(xplot, yplot, color='k', alpha = 0.5, linestyle='dashed')
            else: # gamma 0 (boundaries) highlighted in blue
                ax.plot(xplot, yplot, color='k')

        if plot_frontline:
            for fl_char in self.frontline_characteristics:  # plot frontline points

                ax.plot([fl_char.origin.pos[0]], [fl_char.origin.pos[1]], color='r', linestyle='dashed', marker = 'x')


        # ax.set_ylim(0, 2)
        ax.set_xlim(0, max([p.pos[0] for p in self.frontline_points]))

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        # ax.set_aspect('equal')

        ax.grid()

        if markers: # remove markers if flag is unset
            for line in ax.lines:
                line.set_marker('x')
        if save:
            curdir = os.getcwd()
            plt.savefig(
                "{}/plots/{}_{}.svg".format(curdir, 'geometry', self.iter)
            )
        else:
            plt.show()

        plt.clf()
        plt.close(fig)

    def plot_contours(self, property : str, save = False, plot_characteristics=True, plot_frontline=True, plot_boundaries=True):
        # getattr(a, property, 'default value')

        all_points = self.dead_points.union(self.frontline_points)
        x = np.array(
            [p.pos[0] for p in all_points]
        )

        y = np.array(
            [p.pos[1] for p in all_points]
        )

        z = np.array(
            [getattr(p, property) for p in all_points]
        )  # contourplot value

        fig, ax = plt.subplots(figsize = (8, 6))
        if plot_characteristics or plot_boundaries:
            for dead_char in self.dead_characteristics:

                if dead_char.end is None:
                    continue # if, for some reason, we reach a char that has no end point yet, skip it

                xplot = [dead_char.origin.pos[0], dead_char.end.pos[0]]
                yplot = [dead_char.origin.pos[1], dead_char.end.pos[1]]

                if dead_char.type in [1, -1]: # normal chars colored in grey, other configs can be passed through the advance function
                    if plot_characteristics:
                        ax.plot(xplot, yplot, color='k', alpha = 0.5, linestyle='dashed')
                elif plot_boundaries: # gamma 0 (boundaries) highlighted in blue
                    ax.plot(xplot, yplot, color='k')

        if plot_frontline:
            for fl_char in self.frontline_characteristics: # plot frontline points

                ax.plot([fl_char.origin.pos[0]], [fl_char.origin.pos[1]], color='k', linestyle='dashed', marker = 'x')

        if any(z < 0):
            cmap = 'rainbow'
            max_abs = np.max(np.abs(z))
            vmin = -max_abs
            vmax = max_abs
        else:
            cmap = 'rainbow'
            vmin = np.min(z)
            vmax = np.max(z)

        lvl = sorted(list(set(z)), key=lambda x: x)
        contour = ax.tricontourf(x, y, z, cmap=cmap, levels=lvl, extend='both')
        cbar = fig.colorbar(contour, ax=ax, orientation='vertical', pad=0.1, ticks=np.linspace(vmin, vmax, 10))
        cbar.set_label(property, rotation=90)

        ax.set_xlim(0, max([p.pos[0] for p in self.frontline_points]))

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        # ax.set_aspect('equal')


        ax.grid()

        if save:
            curdir = os.getcwd()
            plt.savefig(
                "{}/plots/{}_{}.svg".format(curdir, property, self.iter)
            )
        else:
            plt.show()

        plt.clf()
        plt.close(fig)




if __name__=="__main__":
    from expansionFan import JetExpansionFan
    from helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
    pm_angle_inlet = prandtl_meyer_from_mach(2.0, 1.4)
    inlet_conditions = GenericFlowElement(pm_angle_inlet, pm_angle_inlet)

    jef = JetExpansionFan(inlet=inlet_conditions, pressure_ratio=1.2, origin=(0, 1), NCHAR=10, gamma=1.4, type=-1)

    inlet_points = [
        FluidPoint((0, yp), pm_angle_inlet, pm_angle_inlet, boundary="lower" if yp==0 else None) for yp in np.linspace(0, 1, 5, endpoint=False)
    ]


    inlet_points.extend(jef.characteristic_origins)

    gc = GeometryCluster(inlet_points)
    gc.run(printFlag=True)
    gc.plot_geometry(save=True)


