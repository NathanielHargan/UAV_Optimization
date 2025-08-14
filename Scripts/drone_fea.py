import numpy as np
import math
from Scripts.beams import BeamSystem
from Scripts.beams import BeamSystem

class DroneFEA:
    def __init__(self, DroneGeometry):
        self.geometry = DroneGeometry
        self.name = DroneGeometry.name
        self.beam_system = BeamSystem(self.name + "_system")
        self.arm_beam = self.geometry.arm_beam
        self.strut_beam = self.geometry.strut_beam
        self.hub_beam = self.geometry.hub_beam

    def create_drone_slice_nodes(self):
        # Origin
        # self.beam_system.add_node(np.array([0, 0, 228.6]), 'center_node')
        self.beam_system.add_node(np.array([0, 0, 0]), 'center_node')
        self.beam_system.add_node(self.geometry.hub_nodes_coords[0], 'hub_node')
        self.beam_system.add_node(self.geometry.strut_nodes_coords[0], 'strut_node')
        self.beam_system.add_node(self.geometry.outer_nodes_coords[0], 'outer_node')
        self.beam_system.add_node(self.geometry.strut_centroid_coords[0], 'strut_node_top')
        self.beam_system.add_node(self.geometry.strut_centroid_coords[-1], 'strut_node_bottom')

    def create_drone_slice_beams(self):

        self.beam_system.add_beam(
            self.hub_beam,
            "center_node",
            "hub_node",
            [0, 0, 1],
            "hub_beam")

        self.beam_system.add_beam(
            self.arm_beam,
            "hub_node",
            "strut_node",
            [0, 0, 1],
            "center_beam")

        self.beam_system.add_beam(
            self.arm_beam,
            "strut_node",
            "outer_node",
            [0, 0, 1],
            "outer_beam")

        self.beam_system.add_beam(
            self.strut_beam,
            "strut_node",
            "strut_node_top",
            [0, 0, 1],
            "strut_element_top")

        self.beam_system.add_beam(
            self.strut_beam,
            "strut_node",
            "strut_node_bottom",
            [0, 0, 1],
            "strut_element_bottom")

    def boundary_conditions_slice(self):
        self.beam_system.add_boundary_condition(0, "en castre", "center_node")

        local_dir_top = self.geometry.strut_centroid_coords[0] - self.geometry.strut_nodes_coords[0]
        local_dir_bot = self.geometry.strut_centroid_coords[-1] - self.geometry.strut_nodes_coords[0]
        self.beam_system.add_boundary_condition(0, "x symm", "strut_node_top", local_dir_top, np.array([0, 0, 1]))
        self.beam_system.add_boundary_condition(0, "x symm", "strut_node_bottom",  local_dir_bot, np.array([0, 0, 1]))

        self.beam_system.add_boundary_condition(0, "y", "outer_node")
        self.beam_system.add_boundary_condition(0, "theta z", "outer_node")
        self.beam_system.add_boundary_condition(0, "theta x", "outer_node")

        self.beam_system.add_boundary_condition(0, "y", "strut_node")
        self.beam_system.add_boundary_condition(0, "theta z", "strut_node")
        self.beam_system.add_boundary_condition(0, "theta x", "strut_node")


        self.beam_system.add_boundary_condition(0, "y", "hub_node")
        self.beam_system.add_boundary_condition(0, "theta z", "hub_node")
        self.beam_system.add_boundary_condition(0, "theta x", "hub_node")


    def create_drone_nodes(self, nominal_rad, strut_pos, blade_num):
        self.nominal_rad = nominal_rad
        self.strut_pos = strut_pos
        self.blade_num = blade_num
        num = 0

        # Origin
        self.beam_system.add_node(
            np.array([0,
            0,
            0]),
            'origin')

        # Rotational Symmetry Nodes
        for theta in np.linspace(0, 2*math.pi*(blade_num-1)/blade_num, num=blade_num):
            # Adds the node where the blades are
            self.beam_system.add_node(
                np.array([math.cos(theta) * nominal_rad,
                math.sin(theta) * nominal_rad,
                0]),
                'outer_node_' + str(num))

            # Adds the node the struts connect to
            self.beam_system.add_node(
                np.array([math.cos(theta) * strut_pos,
                math.sin(theta) * strut_pos,
                0]),
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
        self.beam_system.add_boundary_condition(0, "en castre", "origin")

    def solve_fea(self):
        self.beam_system.solve_FEA()

    def solve_static(self):
        self.beam_system.solve_static()

    def solve_natural_frequencies(self):
        self.beam_system.solve_natural_frequencies()

    def solve_failure(self):
        self.beam_system.solve_failure()

    def solve_failure_swt(self):
        self.beam_system.solve_failure_swt()

    def solve_dynamic_harmonic(self, alpha, beta):
        self.beam_system.solve_dynamic_harmonic(alpha, beta)
