import numpy as np
import math
from Scripts.cross_section_properties import cross_section_circle
from Scripts.cross_section_properties import cross_section_annulus
from Scripts.material_properties import Materials


class BeamSystem:
    def __init__(self, name):
        self.system_name = name  # A string with the name of the system
        self.nodes = np.empty((0, 3))  # List of nodes coordinates
        self.beam_node_indexes = np.empty((0, 2))  # List of the node indexes that beams go between

        self.node_names = []  # List of nodes names
        self.beams = []  # List of beam objects
        self.boundary_conditions_node_index = []  # List of boundary condition objects
        self.boundary_conditions_type = []

        # placeholder
        self.stiffness_matrix = np.empty(0)
        self.f_vector = np.empty(0)

    def add_node(self, x, y, z, name=None):
        coords = np.array([[x, y, z]])
        self.nodes = np.concatenate((self.nodes, coords), axis=0)

        if name is None:
            name = "Node_" + str(len(self.nodes))

        self.node_names.append(name)

    # Function the finds the index of the referred node.
    def select_node(self, node_ref):
        if type(node_ref) is int:
            return node_ref
        elif type(node_ref) is str:
            return self.node_names.index(node_ref)
        else:
            return "error"

    def add_bifurcated_beam(self, beam_type, start_node, end_node, **kwargs):
        if 'start_beam_name' in kwargs:
            start_beam_name = kwargs['start_beam_name']
        else:
            start_beam_name = beam_type.name + "Beam_" + str(len(self.beams))

        if 'end_beam_name' in kwargs:
            end_beam_name = kwargs['end_beam_name']
        else:
            end_beam_name = "Beam_" + str(len(self.beams) + 1)

        if 'midpoint_node_name' in kwargs:
            midpoint_node_name = kwargs['midpoint_node_name']
        else:
            midpoint_node_name = "Node_" + str(len(self.nodes))

        midpoint_node = (self.nodes[self.select_node(start_node)] + self.nodes[self.select_node(end_node)]) / 2
        self.add_node(midpoint_node[0], midpoint_node[1], midpoint_node[2], midpoint_node_name)
        self.add_beam(beam_type, start_node, midpoint_node_name, start_beam_name)
        self.add_beam(beam_type, midpoint_node_name, end_node, end_beam_name)

    def add_beam(self, beam_type, start_node, end_node, name=None):
        start_node_index = self.select_node(start_node)
        end_node_index = self.select_node(end_node)

        beam_node_index = np.array([[start_node_index, end_node_index]])

        self.beam_node_indexes = np.concatenate((self.beam_node_indexes, beam_node_index), axis=0)

        start_node_coords = self.nodes[start_node_index]
        end_node_coords = self.nodes[start_node_index]

        if name is None:
            name = "Beam_" + str(len(self.beams))

        new_beam = Beam(beam_type, start_node_coords, end_node_coords, name)
        self.beams.append(new_beam)

    def create_boundary_condition(self, node_def, bc_type):
        node_index = self.select_node(node_def)
        self.boundary_conditions_node_index.append(node_index)
        self.boundary_conditions_type.append(bc_type)

    def rotate_beam_system(self, rot_x, rot_y, rot_z):
        # Rotation matrix around x-axis
        r_x = np.array([
            [1, 0, 0],
            [0, math.cos(rot_x), -math.sin(rot_x)],
            [0, math.sin(rot_x), math.cos(rot_x)]])

        # Rotation matrix around y-axis
        r_y = np.array([
            [math.cos(rot_y), 0, math.sin(rot_y)],
            [0, 1, 0],
            [-math.sin(rot_y), 0, math.cos(rot_y)]])

        # Rotation matrix around z-axis
        r_z = np.array([
            [math.cos(rot_z),  -math.sin(rot_z), 0],
            [math.sin(rot_z), math.cos(rot_z), 0],
            [0, 0, 1]])

        rotation_matrix = r_x @ r_y @ r_z  # multiplies them all together

        rotated_nodes = np.empty((0, 3))
        for node in self.nodes:
            rotated_node = [rotation_matrix @ np.transpose(node)]
            rotated_nodes = np.concatenate((rotated_nodes, rotated_node), axis=0)

        self.nodes = rotated_nodes  # update nodes

    def create_stiffness_matrix(self):
        stiffness_matrix_size = len(self.nodes) * 2 - len(self.boundary_conditions_node_index)
        self.stiffness_matrix = np.zeros((stiffness_matrix_size,stiffness_matrix_size))
        self.f_vector = np.array(stiffness_matrix_size)


class Beam:
    def __init__(self, beam_type, start_node_coords, end_node_coords, name=None):
        self.name = name  # String with the beam name
        self.beam_type = beam_type

        # Coordinates of the nodes
        self.start_node_coords = np.array(start_node_coords)
        self.end_node_coords = np.array(end_node_coords)

        self.length = np.linalg.norm(np.subtract(end_node_coords, start_node_coords))


class BeamType:
    def __init__(self, cross_section, cross_section_parameters, material, name=None):
        self.cross_section = cross_section
        self.cross_section_parameters = cross_section_parameters
        self.material = Materials[material]

        # Name is optional so
        if name is None:
            self.name = "Beam_Type_" + cross_section + "_" + material
        else:
            self.name = name

        self.cross_section_properties = self.cross_section_properties_init()

    def cross_section_properties_init(self):
        if self.cross_section.lower() == "circle":
            return cross_section_circle(self.cross_section_parameters)
        elif self.cross_section.lower() == "annulus":
            return cross_section_annulus(self.cross_section_parameters)
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'





