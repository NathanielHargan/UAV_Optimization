import numpy as np
import math
from beams import BeamSystem


class Drone:
    def __init__(self, name):
        self.name = name # Name of the drone
        self.beam_system = BeamSystem(self.name + "_system")

    def create_drone_beams(self, hub_rad, nominal_rad, strut_pos, blade_num):
        num = 0
        # Rotational Symmetry
        for theta in np.linspace(0, 2*math.pi*(blade_num-1)/blade_num, num=blade_num):
            # Adds the node where the blades are
            self.beam_system.add_node(
                math.cos(theta) * nominal_rad,
                math.sin(theta) * nominal_rad,
                0,
                'outer_node_' + str(num))

            # Adds the node the struts connect to
            self.beam_system.add_node(
                math.cos(theta) * strut_pos,
                math.sin(theta) * strut_pos,
                0,
                "strut_node_" + str(num))

            # Adds the node that connects the arms to the hub
            self.beam_system.add_node(
                math.cos(theta) * hub_rad,
                math.sin(theta) * hub_rad,
                0,
                "hub_node_" + str(num))

            num += 1

        self.beam_system.add_beam(
            1,
            2,
            "annulus",
            [0.1,0.2])














    #%%
