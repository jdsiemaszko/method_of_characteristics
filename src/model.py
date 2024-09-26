from fluidPoint import FluidPoint
from boundary import BoundaryLine
import numpy as np
from helper import prandtl_meyer_from_mach, mach_from_prandtl_meyer
import matplotlib.pyplot as plt
from characteristic import Characteristic
from expansion_fan import JetExpansionFan

class Model():
    def __init__(self, config):

        self.config = config

        self.characteristics = [] # initialize characteristics as empty

    def run(self):

        self.expansion_fan(self.point_A)

    def create_inlet_characteristics(self):
        for yposition in




