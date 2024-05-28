import numpy as np
import math
import scipy
from Scripts.cross_section_properties import cross_section_circle
from Scripts.cross_section_properties import cross_section_annulus
from Scripts.material_properties import Materials
from Scripts.material_properties import Gravity
import Scripts.FEA_3D as FEA_3D
import time
import matplotlib.pyplot as plt


class BeamSystem:
    def __init__(self, name):
        self.system_name = name  # A string with the name of the system
        self.nodes = np.empty((0, 3))  # List of nodes coordinates
        self.beam_node_indexes = np.empty((0, 2))  # List of the node indexes that beams go between

        self.forces = np.empty((0, 3))  # List of Vectors
        self.force_node_indexes = np.empty((0, 1)) # List of node indexes

        # The value of the BC
        self.boundary_conditions = np.empty((0, 1))
        # x | y | z | theta x | theta y | theta z  alternatively 0 | 1 | 2 | 3 | 4 | 5
        self.boundary_conditions_type = np.empty((0, 1))
        # The boundary condition node index vector
        self.boundary_conditions_node_indexes = np.empty((0, 1))

        self.node_names = []  # List of nodes names
        self.beams = []  # List of beam objects

        # placeholder
        self.global_stiffness_matrix = np.empty(0)
        self.ru = np.empty(0)
        self.dp = np.empty(0)
        self.p_index = np.empty(0)
        self.u_index = np.empty(0)
        self.displacement_angle_vector = np.empty(0)
        self.force_moment_vector = np.empty(0)
        self.kuu = np.empty(0)
        self.kup = np.empty(0)
        self.kpu = np.empty(0)
        self.kpp = np.empty(0)
        self.x_displacements = np.empty(0)
        self.y_displacements = np.empty(0)
        self.z_displacements = np.empty(0)
        self.x_angles = np.empty(0)
        self.y_angles = np.empty(0)
        self.z_angles = np.empty(0)
        self.x_forces = np.empty(0)
        self.y_forces = np.empty(0)
        self.z_forces = np.empty(0)
        self.x_moments = np.empty(0)
        self.y_moments = np.empty(0)
        self.z_moments = np.empty(0)

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

    def select_boundary_condition_type(self, bc_ref):
        if type(bc_ref) is int:
            return bc_ref
        elif type(bc_ref) is str:
            return ["x", "y", "z", "theta x", "theta y", "theta z"].index(bc_ref)  # life hack
        else:
            return "error"

    def add_boundary_condition(self, boundary_val, bc_type_ref, node_ref):
        boundary_val_ = np.array([[boundary_val]])
        node_index = np.array([[self.select_node(node_ref)]])
        bc_type = np.array([[self.select_boundary_condition_type(bc_type_ref)]])
        self.boundary_conditions = np.concatenate((self.boundary_conditions, boundary_val_), axis=0)
        self.boundary_conditions_type = np.concatenate((self.boundary_conditions_type, bc_type), axis=0)
        self.boundary_conditions_node_indexes = np.concatenate((self.boundary_conditions_node_indexes, node_index), axis=0)

    def add_force(self, x, y, z, node_ref):
        node = np.array([[self.select_node(node_ref)]])
        vector = np.array([[x, y, z]])
        self.forces = np.concatenate((self.forces, vector), axis=0)
        self.force_node_indexes = np.concatenate((self.force_node_indexes, node), axis=0)
    '''
    def add_bifurcated_beam(self, beam_type, start_node, end_node, k_node, **kwargs):
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
        self.add_beam(beam_type, start_node, midpoint_node_name, k_node, start_beam_name)
        self.add_beam(beam_type, midpoint_node_name, end_node, k_node, end_beam_name)

    def add_point_mass_beam(self, beam_type, start_node, end_node, k_node, **kwargs):
        self.add_bifurcated_beam(beam_type, start_node, end_node, k_node, **kwargs)
        length = np.linalg.norm(np.subtract(
            self.nodes[self.select_node(end_node)],
            self.nodes[self.select_node(start_node)]))
        mass = -beam_type.material_properties["density"] * beam_type.cross_section_properties["area"] * length * Gravity
        self.add_force(0, 0, mass, int(len(self.nodes)-1))  # add force on new node
    '''
    def add_beam(self, beam_type, start_node, end_node, k_node, name=None):
        start_node_index = self.select_node(start_node)
        end_node_index = self.select_node(end_node)
        k_node = np.array(k_node)

        beam_node_index = np.array([[start_node_index, end_node_index]])

        self.beam_node_indexes = np.concatenate((self.beam_node_indexes, beam_node_index), axis=0)

        start_node_coords = self.nodes[start_node_index]
        end_node_coords = self.nodes[end_node_index]

        if name is None:
            name = "Beam_" + str(len(self.beams))

        new_beam = Beam(beam_type, start_node_coords, end_node_coords, k_node, name)
        self.beams.append(new_beam)

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

    def solve_FEA(self):
        global_length = np.shape(self.nodes)[0] * 6

        dp = self.boundary_conditions[:, 0]
        p_index = np.intc(self.boundary_conditions_node_indexes * 6 + self.boundary_conditions_type)

        u_index_predel = np.intc(np.arange(global_length))
        ru_predel = np.zeros(global_length)  # fu before deleting bcs

        for i in range(len(self.forces)):  # for loop + counter
            r_index = int(self.force_node_indexes[i][0] * 6)
            ru_predel[r_index] += self.forces[i][0]
            ru_predel[r_index+1] += self.forces[i][1]
            ru_predel[r_index+2] += self.forces[i][2]

        ru = np.delete(ru_predel, p_index, axis=0)
        u_index = np.intc(np.delete(u_index_predel, p_index, axis=0))

        self.dp = dp
        self.p_index = p_index
        self.ru = ru
        self.u_index = u_index

        up_index = np.concatenate((u_index, p_index), axis=None)
        global_element_matrices = np.empty((len(self.beam_node_indexes), 12, 12))

        for i in range(len(self.beam_node_indexes)):
            global_element_matrices[i] = self.beams[i].global_stiffness_matrix

        self.global_stiffness_matrix = FEA_3D.assemble_stiffness_3d(self.beam_node_indexes, global_element_matrices, global_length, up_index)

        kuu, kup, kpu, kpp = FEA_3D.partition_stiffness_matrix(self.global_stiffness_matrix, len(u_index))

        self.kuu = kuu
        self.kup = kup
        self.kpu = kpu
        self.kpp = kpp

        du = scipy.linalg.solve(kuu, ru - (kup @ dp))

        rp = (kpu @ du) + (kpp @ dp)

        dup = np.concatenate((du, dp), axis=None)
        rup = np.concatenate((ru, rp), axis=None)

        self.du = du
        self.rp = rp
        self.rup = rup
        self.dup = dup

        d = np.empty(global_length)
        r = np.empty(global_length)

        for i in range(global_length):
            d[up_index[i]] = dup[i]
            r[up_index[i]] = rup[i]

        self.displacement_angle_vector = d
        self.force_moment_vector = r

        self.x_displacements = d[0::6]
        self.y_displacements = d[1::6]
        self.z_displacements = d[2::6]
        self.x_angles = d[3::6]
        self.y_angles = d[4::6]
        self.z_angles = d[5::6]

        self.x_forces = r[0::6]
        self.y_forces = r[1::6]
        self.z_forces = r[2::6]
        self.x_moments = r[3::6]
        self.y_moments = r[4::6]
        self.z_moments = r[5::6]

        for beam_num, beam in enumerate(self.beams):
            beam_node_0 = int(self.beam_node_indexes[beam_num, 0])
            beam_node_1 = int(self.beam_node_indexes[beam_num, 1])
            beam.x_displacement = np.array([self.x_displacements[beam_node_0], self.x_displacements[beam_node_1]])
            beam.y_displacement = np.array([self.y_displacements[beam_node_0], self.y_displacements[beam_node_1]])
            beam.z_displacement = np.array([self.z_displacements[beam_node_0], self.z_displacements[beam_node_1]])
            beam.x_angles = np.array([self.x_angles[beam_node_0], self.x_angles[beam_node_1]])
            beam.y_angles = np.array([self.y_angles[beam_node_0], self.y_angles[beam_node_1]])
            beam.z_angles = np.array([self.z_angles[beam_node_0], self.z_angles[beam_node_1]])


class Beam:
    def __init__(self, beam_type, start_node_coords, end_node_coords, k_node, name=None):
        self.name = name  # String with the beam name
        self.k_node = k_node
        self.beam_type = beam_type

        # Coordinates of the nodes
        self.start_node_coords = np.array(start_node_coords)
        self.end_node_coords = np.array(end_node_coords)

        self.length = np.linalg.norm(np.subtract(end_node_coords, start_node_coords))

        # a, e, l, g, i_y, i_z, k, k_y, k_z
        self.local_stiffness_matrix = FEA_3D.local_stiffness_3d(
            beam_type.cross_section_properties["area"],
            beam_type.material_properties["elastic modulus"],
            self.length,
            beam_type.material_properties["shear modulus"],
            beam_type.cross_section_properties['second moment of area x'],
            beam_type.cross_section_properties['second moment of area y'],
            beam_type.cross_section_properties["torsional constant"],
            beam_type.cross_section_properties['transverse shear deflection constant x'],
            beam_type.cross_section_properties['transverse shear deflection constant y']

        )

        self.mass_matrix = FEA_3D.mass_matrix_3d(
            beam_type.cross_section_properties["area"],
            self.length,
            beam_type.material_properties["mass density"]
        )

        direction = self.end_node_coords - self.start_node_coords
        self.transformation_matrix = FEA_3D.transformation_matrix(direction, k_node)

        self.global_stiffness_matrix = FEA_3D.local_to_global_stiffness_matrix(
            self.local_stiffness_matrix,
            self.transformation_matrix)

        self.x_displacement = []
        self.y_displacement = []
        self.z_displacement = []
        self.x_angles = []
        self.y_angles = []
        self.z_angles = []


class BeamType:
    def __init__(self, cross_section, cross_section_parameters, material, name=None):
        self.cross_section = cross_section
        self.cross_section_parameters = cross_section_parameters
        self.material_properties = Materials[material]

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





