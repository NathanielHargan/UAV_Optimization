import numpy as np
import math
from Scripts.Finite_Elements.Beam_FEA_System.beam_system import BeamSystem
from Scripts.Finite_Elements.Beam_FEA_System.beam_type import BeamType


class DroneFEA:
    def __init__(self, DroneGeometry):
        self.geometry = DroneGeometry
        self.name = DroneGeometry.name
        self.beam_system = BeamSystem(self.name + "_system")
        self.arm_beam = self.geometry.arm_beam
        self.strut_beam = self.geometry.strut_beam
        self.hub_beam = self.geometry.hub_beam
        self.hub_sections = self.geometry.hub_sections
        self.final_hub_node_name = ""

    def create_drone_slice_nodes(self):
        # Origin
        # self.beam_system.add_node(np.array([0, 0, 228.6]), 'center_node')

        # Add hub nodes
        for i in range(self.hub_sections+1):
            hub_node_name = 'hub_node_' + str(i)
            self.beam_system.add_node(self.geometry.hub_nodes_coords[0, i], hub_node_name)

        self.beam_system.add_node(self.geometry.strut_nodes_coords[0], 'strut_node')
        self.beam_system.add_node(self.geometry.outer_nodes_coords[0], 'outer_node')
        self.beam_system.add_node(self.geometry.strut_centroid_coords[0], 'strut_node_top')
        self.beam_system.add_node(self.geometry.strut_centroid_coords[-1], 'strut_node_bottom')

    def create_drone_slice_beams(self):
        # Add hub beams
        for i in range(self.hub_sections):
            parameters = self.hub_beam.cross_section_parameters
            parameters[0] = self.geometry.beam_widths[i]
            hub_beam_type = BeamType(self.hub_beam.cross_section,parameters,self.hub_beam.material,self.hub_beam.name)
            hub_beam_name = 'hub_beam_' + str(i)
            hub_node_name_0 = 'hub_node_' + str(i)
            hub_node_name_1 = 'hub_node_' + str(i+1)
            self.beam_system.add_beam(
                hub_beam_type,
                hub_node_name_0,
                hub_node_name_1,
                [0, 0, 1],
                hub_beam_name)

        self.final_hub_node_name = 'hub_node_' + str(self.hub_sections)
        self.beam_system.add_beam(
            self.arm_beam,
             self.final_hub_node_name,
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
        self.beam_system.add_boundary_condition(0, "en castre", "hub_node_0")

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


        self.beam_system.add_boundary_condition(0, "y", self.final_hub_node_name)
        self.beam_system.add_boundary_condition(0, "theta z", self.final_hub_node_name)
        self.beam_system.add_boundary_condition(0, "theta x", self.final_hub_node_name)

    def create_drone_nodes(self):
        # Origin
        self.beam_system.add_node(
            np.array([0,
            0,
            0]),
            'origin')
        for i in range(self.geometry.blade_num):
            for j in range(1, self.hub_sections+1):
                hub_node_name = 'hub_node_' + str(j) + '_' + str(i)
                self.beam_system.add_node(self.geometry.hub_nodes_coords[i, j], hub_node_name)

            self.beam_system.add_node(self.geometry.strut_nodes_coords[i], f'strut_node_{i}')
            self.beam_system.add_node(self.geometry.outer_nodes_coords[i], f'outer_node_{i}')

    def create_drone_beams(self):
        arm_beam = self.geometry.arm_beam
        strut_beam = self.geometry.strut_beam

        # Adding Beams
        for i in range(0, self.geometry.blade_num):
            # Outer to Strut
            self.beam_system.add_beam(
                arm_beam,
                f'strut_node_{i}',
                f'outer_node_{i}',
                [0, 0, 1],
                f'outer_arm_beam_{i}')

            # Strut to Hub
            self.beam_system.add_beam(
                arm_beam,
                "strut_node_" + str(i),
                'hub_node_' + str(self.hub_sections) + '_' + str(i),
                [0, 0, 1],
                f'inner_arm_beam_{i}')

            parameters = self.hub_beam.cross_section_parameters
            parameters[0] = self.geometry.beam_widths[0]
            hub_beam_type = BeamType(self.hub_beam.cross_section, parameters, self.hub_beam.material,
                                     self.hub_beam.name)

            self.beam_system.add_beam(
                hub_beam_type,
                "origin",
                'hub_node_1_' + str(i),
                [0, 0, 1],
            f'hub_beam_0_{i}')

            for j in range(self.hub_sections-1):
                parameters = self.hub_beam.cross_section_parameters
                parameters[0] = self.geometry.beam_widths[j+1]
                hub_beam_type = BeamType(self.hub_beam.cross_section, parameters, self.hub_beam.material,
                                         self.hub_beam.name)
                self.beam_system.add_beam(
                    hub_beam_type,
                    'hub_node_' + str(j+1) + '_' + str(i),
                    'hub_node_' + str(j+2) + '_' + str(i),
                    [0, 0, 1],
                f'hub_beam_{j+1}_{i}')


        # Final to 0
        self.beam_system.add_beam(
            strut_beam,
            "strut_node_" + str(self.geometry.blade_num-1),
            "strut_node_0",
            [0, 0, 1],
            f'strut_beam_{self.geometry.blade_num-1}')

        for i in range(0,  self.geometry.blade_num-1):
            # Strut to Strut
            self.beam_system.add_beam(
                strut_beam,
                f'strut_node_{i}',
                f'strut_node_{i+1}',
                [0, 0, 1],
            f'strut_beam_{i}')

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
