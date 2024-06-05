import numpy as np
import math
from Scripts.beams import BeamSystem


class Drone:
    def __init__(self, name):
        self.name = name  # Name of the drone
        self.beam_system = BeamSystem(self.name + "_system")

    def create_drone_slice_nodes(self, nominal_rad, strut_pos, blade_num):
        # Origin
        self.beam_system.add_node(
            0,
            0,
            0,
            'center_node')

        self.beam_system.add_node(
            strut_pos,
            0,
            0,
            'strut_node')

        self.beam_system.add_node(
            nominal_rad,
            0,
            0,
            'outer_node')

        theta = 2 * math.pi / blade_num
        x_next = math.cos(theta) * strut_pos
        y_next = math.sin(theta) * strut_pos
        x = (x_next + strut_pos)/2
        y = y_next/2
        self.beam_system.add_node(
            x,
            y,
            0,
            'strut_node_top')

        self.beam_system.add_node(
            x,
            -y,
            0,
            'strut_node_bottom')

    def create_drone_slice_beams(self, arm_beam, strut_beam):
        self.beam_system.add_beam(
            arm_beam,
            "center_node",
            "strut_node",
            [0, 0, 1])

        self.beam_system.add_beam(
            arm_beam,
            "strut_node",
            "outer_node",
            [0, 0, 1])

        self.beam_system.add_beam(
            strut_beam,
            "strut_node",
            "strut_node_top",
            [0, 0, 1])

        self.beam_system.add_beam(
            strut_beam,
            "strut_node",
            "strut_node_bottom",
            [0, 0, 1])

    def boundary_conditions_slice(self):
        self.beam_system.add_boundary_condition(0, "x", "center_node")
        self.beam_system.add_boundary_condition(0, "y", "center_node")
        self.beam_system.add_boundary_condition(0, "z", "center_node")
        self.beam_system.add_boundary_condition(0, "theta x", "center_node")
        self.beam_system.add_boundary_condition(0, "theta y", "center_node")
        self.beam_system.add_boundary_condition(0, "theta z", "center_node")

        self.beam_system.add_boundary_condition(0, "x", "strut_node_top")
        self.beam_system.add_boundary_condition(0, "y", "strut_node_top")

        self.beam_system.add_boundary_condition(0, "x", "strut_node_bottom")
        self.beam_system.add_boundary_condition(0, "y", "strut_node_bottom")



    def create_drone_nodes(self, nominal_rad, strut_pos, blade_num):
        self.nominal_rad = nominal_rad
        self.strut_pos = strut_pos
        self.blade_num = blade_num
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

    def create_drone_beams(self, arm_beam, strut_beam):
        # Adding Beams
        for i in range(0, self.blade_num):
            # Outer to Strut
            self.beam_system.add_beam(
                arm_beam,
                "outer_node_" + str(i),
                "strut_node_" + str(i),
                [0, 0, 1])

            # Strut to Hub
            self.beam_system.add_beam(
                arm_beam,
                "strut_node_" + str(i),
                "origin",
                [0, 0, 1])

        # Final to 0
        self.beam_system.add_beam(
            strut_beam,
            "strut_node_" + str(self.blade_num-1),
            "strut_node_0",
            [0, 0, 1])

        for i in range(0,  self.blade_num-1):
            # Strut to Strut
            self.beam_system.add_beam(
                strut_beam,
                "strut_node_" + str(i),
                "strut_node_" + str(i+1),
                [0, 0, 1])

    def rotate_drone_nodes(self, x, y, z):
        self.beam_system.rotate_beam_system(x, y, z)

    def boundary_conditions(self):
        self.beam_system.add_boundary_condition(0, "x", "origin")
        self.beam_system.add_boundary_condition(0, "y", "origin")
        self.beam_system.add_boundary_condition(0, "z", "origin")
        self.beam_system.add_boundary_condition(0, "theta x", "origin")
        self.beam_system.add_boundary_condition(0, "theta y", "origin")
        self.beam_system.add_boundary_condition(0, "theta z", "origin")

    def solve_fea(self):
        self.beam_system.solve_FEA()













    #%%

#%%
