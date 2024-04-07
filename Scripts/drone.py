import numpy as np
import math
from Scripts.beams import BeamSystem


class Drone:
    def __init__(self, name):
        self.name = name  # Name of the drone
        self.beam_system = BeamSystem(self.name + "_system")

    def create_drone(self, nominal_rad, strut_pos, blade_num, arm_beam, strut_beam):

        num = 0

        # Origin
        self.beam_system.add_node(
            0,
            0,
            0,
            'origin')

        # Rotational Symmetry Nodes
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
            num += 1

        # Adding Beams
        for i in range(0, blade_num):
            # Outer to Strut
            self.beam_system.add_beam(
                arm_beam,
                "outer_node_" + str(i),
                "strut_node_" + str(i))

            # Strut to Hub
            self.beam_system.add_beam(
                arm_beam,
                "strut_node_" + str(i),
                "origin")

        # Final to 0
        self.beam_system.add_beam(
            strut_beam,
            "strut_node_" + str(blade_num-1),
            "strut_node_0")

        for i in range(0, blade_num-1):
            # Strut to Strut
            self.beam_system.add_beam(
                strut_beam,
                "strut_node_" + str(i),
                "strut_node_" + str(i+1))















    #%%
