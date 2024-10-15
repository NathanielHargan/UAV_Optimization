#%%
import numpy as np
import math
import scipy
from Scripts.cross_section_properties import cross_section_circle
from Scripts.cross_section_properties import cross_section_annulus
from Scripts.cross_section_properties import cross_section_rectangle
from Scripts.cross_section_properties import cross_section_hexagon
from Scripts.material_properties import Materials
from Scripts.material_properties import Gravity
import Scripts.FEA_3D as FEA_3D
import time
import matplotlib.pyplot as plt

import logging
from numpy.testing import assert_almost_equal


class Node:
    def __init__(self, location, node_index, name):
        self.location = location
        self.node_index = node_index
        self.name = name

        self.force = np.array([0, 0, 0])
        self.moment = np.array([0, 0, 0])
        self.local_bc_transform = np.identity(6)
        self.result_displacement = np.array([0, 0, 0])
        self.result_angular_displacement = np.array([0, 0, 0])
        self.result_force = np.array([0, 0, 0])
        self.result_moment = np.array([0, 0, 0])

    def add_force(self, force):
        self.force = self.force + force

    def add_moment(self, moment):
        self.moment = self.moment + moment
        # To do:
        # Make node class
        # Move more information to the beam class
        # Local transformation option on elements for a specified node.

    def create_local_transform(self, transform_dir, k_node_dir):
        self.local_bc_transform = FEA_3D.transformation_matrix_node(transform_dir, k_node_dir)


class Beam:
    def __init__(self, beam_type, start_node, end_node, k_node, name=None):
        self.name = name  # String with the beam name
        self.k_node = k_node
        self.beam_type = beam_type
        self.start_node = start_node
        self.end_node = end_node

        self.length = np.linalg.norm(np.subtract(end_node.location, start_node.location))
        self.direction = end_node.location - start_node.location

        self.d_local = np.empty(12)
        self.r_local = np.empty(12)

        self.local_stiffness_matrix = np.empty((12, 12))
        self.mass_matrix = np.empty((12, 12))
        self.local_bc_transformation = np.identity(12)
        self.transformation_matrix = np.empty((12, 12))
        self.global_stiffness_matrix = np.empty((12, 12))
        self.applied_stiffness_matrix = np.empty((12, 12))

        self.x_displacements_local = np.empty(2)
        self.y_displacements_local = np.empty(2)
        self.z_displacements_local = np.empty(2)
        self.x_angles_local = np.empty(2)
        self.y_angles_local = np.empty(2)
        self.z_angles_local = np.empty(2)
        self.x_forces_local = np.empty(2)
        self.y_forces_local = np.empty(2)
        self.z_forces_local = np.empty(2)
        self.x_moments_local = np.empty(2)
        self.y_moments_local = np.empty(2)
        self.z_moments_local = np.empty(2)

        self.x_displacements_global = np.empty(2)
        self.y_displacements_global = np.empty(2)
        self.z_displacements_global = np.empty(2)
        self.x_angles_global = np.empty(2)
        self.y_angles_global = np.empty(2)
        self.z_angles_global = np.empty(2)
        self.x_forces_global = np.empty(2)
        self.y_forces_global = np.empty(2)
        self.z_forces_global = np.empty(2)
        self.x_moments_global = np.empty(2)
        self.y_moments_global = np.empty(2)
        self.z_moments_global = np.empty(2)

        self.a_inv_y_shape_vector = np.empty(4)
        self.a_inv_z_shape_vector = np.empty(4)

        self.solve_stiffness_matrix()

        self.stresses_bending_y = np.empty(2)
        self.stresses_bending_z = np.empty(2)
        self.stresses_axial = np.empty(2)

    def solve_stiffness_matrix(self):
        # a, e, l, g, i_y, i_z, k, k_y, k_z
        self.local_stiffness_matrix = FEA_3D.local_stiffness_3d(
            self.beam_type.cross_section_properties["area"],
            self.beam_type.material_properties["elastic modulus"],
            self.length,
            self.beam_type.material_properties["shear modulus"],
            self.beam_type.cross_section_properties['second moment of area y'],
            self.beam_type.cross_section_properties['second moment of area z'],
            self.beam_type.cross_section_properties["torsional constant"],
            self.beam_type.cross_section_properties['transverse shear deflection constant y'],
            self.beam_type.cross_section_properties['transverse shear deflection constant z']
        )

        self.mass_matrix = FEA_3D.mass_matrix_3d(
            self.beam_type.cross_section_properties["area"],
            self.length,
            self.beam_type.material_properties["mass density"]
        )

        self.transformation_matrix = FEA_3D.transformation_matrix_element(self.direction, self.k_node)

        self.global_stiffness_matrix = FEA_3D.local_to_global_stiffness_matrix(
            self.local_stiffness_matrix,
            self.transformation_matrix) #  np.transpose(t) @ k @ t


    def solve_local(self, d, r):
        self.d_local = self.transformation_matrix @ d
        self.r_local = self.transformation_matrix @ r
        self.x_displacements_local = np.array([self.d_local[0], self.d_local[6]])
        self.y_displacements_local = np.array([self.d_local[1], self.d_local[7]])
        self.z_displacements_local = np.array([self.d_local[2], self.d_local[8]])
        self.x_angles_local = np.array([self.d_local[3], self.d_local[9]])
        self.y_angles_local = np.array([self.d_local[4], self.d_local[10]])
        self.z_angles_local = np.array([self.d_local[5], self.d_local[11]])

        self.x_forces_local = np.array([self.r_local[0], self.r_local[6]])
        self.y_forces_local = np.array([self.r_local[1], self.r_local[7]])
        self.z_forces_local = np.array([self.r_local[2], self.r_local[8]])
        self.x_moments_local = np.array([self.r_local[3], self.r_local[9]])
        self.y_moments_local = np.array([self.r_local[4], self.r_local[10]])
        self.z_moments_local = np.array([self.r_local[5], self.r_local[11]])

        self.stresses_bending_y = self.y_moments_local /self.beam_type.cross_section_properties['second moment of area y']
        self.stresses_bending_z = self.z_moments_local /self.beam_type.cross_section_properties['second moment of area z']
        self.stresses_axial = self.y_forces_local / self.beam_type.cross_section_properties["area"]


    def solve_shape_functions(self):
        a_inv_y = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_y)
        a_inv_z = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_z)

        self.a_inv_y_shape_vector = a_inv_y @ np.array([self.y_displacements_local[0],
                                                        self.z_angles_local[0],
                                                        self.y_displacements_local[1],
                                                        self.z_angles_local[1]])

        self.a_inv_z_shape_vector = a_inv_z @ np.array([self.z_displacements_local[0],
                                                        self.y_angles_local[0],
                                                        self.z_displacements_local[1],
                                                        self.y_angles_local[1]])

    def return_shape_functions(self, xi):
        x = xi * self.length
        X_y = np.array([
            [1, x, x**2, x**3],
            [0, 1, 2*x, (3 * x ** 2) - 6 * self.beam_type.g_y],
            [0, 0, 0, -6 * self.beam_type.g_y]])

        X_z = np.array([
            [1, x, x**2, x**3],
            [0, 1, 2*x, (3 * x ** 2) - 6 * self.beam_type.g_z],
            [0, 0, 0, -6 * self.beam_type.g_z]])

        y = X_y @ self.a_inv_y_shape_vector
        z = X_z @ self.a_inv_z_shape_vector

        return y, z


class BeamSystem:
    def __init__(self, name):
        self.system_name = name  # A string with the name of the system
        self.beam_node_indexes = np.empty((0, 2))  # List of the node indexes that beams go between

        self.boundary_conditions_node_indexes = np.empty((0, 1))
        self.boundary_conditions_type = np.empty((0, 1))
        self.boundary_conditions = np.empty((0, 1))

        self.nodes = []  # list of node objects
        self.beams = []  # List of beam objects

        # list of node names for calling
        self.node_names = []  # list of node names
        self.beam_names = []  # list of beam names

        # placeholder
        self.global_stiffness_matrix = np.empty(0)
        self.local_bc_transform = np.empty(0)
        self.transformed_stiffness_matrix = np.empty(0)
        self.applied_stiffness_matrix = np.empty(0)

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


    def add_node(self, location, name=None):
        if name is None:
            name = "Node_" + str(len(self.nodes))

        node_index = len(self.nodes)
        self.nodes.append(Node(location, node_index, name))
        self.node_names.append(name)

    # Function the finds the index of the referred node.
    def select_node(self, node_ref):
        if type(node_ref) is int:
            return node_ref
        elif type(node_ref) is str:
            return self.node_names.index(node_ref)
        else:
            return "error"

    def select_element(self, element_ref):
        if type(element_ref) is int:
            return element_ref
        elif type(element_ref) is str:
            return self.beam_names.index(element_ref)
        else:
            return "error"

    def select_boundary_condition_type(self, bc_ref):
        if type(bc_ref) is int:
            return bc_ref
        elif type(bc_ref) is str:
            return ["x", "y", "z", "theta x", "theta y", "theta z"].index(bc_ref)  # life hack
        else:
            print("invalid boundary condition type")
            return "error"

    def add_boundary_condition(self, boundary_val, bc_type_ref, node_ref, direction=None, k_node_dir=None):
        node_index = self.select_node(node_ref)
        bc_type = self.select_boundary_condition_type(bc_type_ref)
        self.boundary_conditions = np.concatenate((self.boundary_conditions, [[boundary_val]]))
        self.boundary_conditions_type = np.concatenate((self.boundary_conditions_type, [[bc_type]]))
        self.boundary_conditions_node_indexes = np.concatenate((self.boundary_conditions_node_indexes, [[node_index]]))

        if direction is not None:
            node = self.nodes[node_index]
            node.create_local_transform(direction, k_node_dir)

    def add_force(self, direction, node_ref):
        node = self.nodes[self.select_node(node_ref)]
        node.add_force(direction)

    def add_moment(self, direction, node_ref):
        node = self.nodes[self.select_node(node_ref)]
        node.add_moment(direction)

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

        if name is None:
            name = "Beam_" + str(len(self.beams))

        new_beam = Beam(beam_type, self.nodes[start_node_index], self.nodes[end_node_index], k_node, name)
        self.beams.append(new_beam)
        self.beam_names.append(name)
    '''
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
    '''

    def solve_FEA(self):
        global_length = np.shape(self.nodes)[0] * 6

        dp = self.boundary_conditions[:, 0]
        p_index = np.intc(self.boundary_conditions_node_indexes * 6 + self.boundary_conditions_type)

        u_index_predel = np.intc(np.arange(global_length))
        ru_predel = np.zeros(global_length)  # fu before deleting bcs

        for i in range(len(self.nodes)):  # for loop + counter
            r_index = int(i * 6)
            ru_predel[r_index] += self.nodes[i].force[0]
            ru_predel[r_index+1] += self.nodes[i].force[1]
            ru_predel[r_index+2] += self.nodes[i].force[2]

        for i in range(len(self.nodes)):  # for loop + counter
            r_index = int(i * 6)
            ru_predel[r_index+3] += self.nodes[i].moment[0]
            ru_predel[r_index+4] += self.nodes[i].moment[1]
            ru_predel[r_index+5] += self.nodes[i].moment[2]

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

        self.global_stiffness_matrix = FEA_3D.assemble_stiffness_3d(self.beam_node_indexes, global_element_matrices, global_length)

        local_bc_transform = np.zeros([global_length, global_length])

        for i, node in enumerate(self.nodes):
            local_bc_transform[i*6:6+i*6, i*6:6+i*6] = node.local_bc_transform

        self.local_bc_transform = local_bc_transform

        self.transformed_stiffness_matrix = self.local_bc_transform @ self.global_stiffness_matrix @ np.transpose(self.local_bc_transform)

        self.applied_stiffness_matrix = FEA_3D.rearrange_stiffness_matrix(self.transformed_stiffness_matrix, up_index)

        kuu, kup, kpu, kpp = FEA_3D.partition_stiffness_matrix(self.applied_stiffness_matrix, len(u_index))

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

        self.displacement_angle_vector = np.transpose(self.local_bc_transform) @ d
        self.force_moment_vector = np.transpose(self.local_bc_transform) @ r

        self.x_displacements = self.displacement_angle_vector[0::6]
        self.y_displacements = self.displacement_angle_vector[1::6]
        self.z_displacements = self.displacement_angle_vector[2::6]

        self.x_angles = self.displacement_angle_vector[3::6]
        self.y_angles = self.displacement_angle_vector[4::6]
        self.z_angles = self.displacement_angle_vector[5::6]

        self.x_forces = self.force_moment_vector[0::6]
        self.y_forces = self.force_moment_vector[1::6]
        self.z_forces = self.force_moment_vector[2::6]

        self.x_moments = self.force_moment_vector[3::6]
        self.y_moments = self.force_moment_vector[4::6]
        self.z_moments = self.force_moment_vector[5::6]

        for i, node in enumerate(self.nodes):
            self.nodes[i].result_displacement = np.array([
                self.x_displacements[i],
                self.y_displacements[i],
                self.z_displacements[i]])

            self.nodes[i].result_angular_displacement = np.array([
                self.x_angles[i],
                self.y_angles[i],
                self.z_angles[i]])

            self.nodes[i].result_force = np.array([
                self.x_forces[i],
                self.y_forces[i],
                self.z_forces[i]])

            self.nodes[i].result_moment = np.array([
                self.x_moments[i],
                self.y_moments[i],
                self.z_moments[i]])

        for beam_num, beam in enumerate(self.beams):
            beam_index_0 = int(self.beam_node_indexes[beam_num][0] * 6)
            beam_index_1 = int(self.beam_node_indexes[beam_num][1] * 6)
            local_d = np.concatenate((d[beam_index_0:beam_index_0 + 6], d[beam_index_1:beam_index_1+6]))
            local_r = np.concatenate((r[beam_index_0:beam_index_0 + 6], r[beam_index_1:beam_index_1+6]))
            beam.solve_local(local_d, local_r)

            beam_node_index_0 = int(self.beam_node_indexes[beam_num][0])
            beam_node_index_1 = int(self.beam_node_indexes[beam_num][1])

            beam.x_displacements_global = np.array([
                self.x_displacements[beam_node_index_0],
                self.x_displacements[beam_node_index_1]])
            beam.y_displacements_global = np.array([
                self.y_displacements[beam_node_index_0],
                self.y_displacements[beam_node_index_1]])
            beam.z_displacements_global = np.array([
                self.z_displacements[beam_node_index_0],
                self.z_displacements[beam_node_index_1]])

            beam.x_angles_global = np.array([
                self.x_angles[beam_node_index_0],
                self.x_angles[beam_node_index_1]])
            beam.y_angles_global = np.array([
                self.y_angles[beam_node_index_0],
                self.y_angles[beam_node_index_1]])
            beam.z_angles_global = np.array([
                self.z_angles[beam_node_index_0],
                self.z_angles[beam_node_index_1]])

            beam.x_forces_global = np.array([
                self.x_forces[beam_node_index_0],
                self.x_forces[beam_node_index_1]])
            beam.y_forces_global = np.array([
                self.y_forces[beam_node_index_0],
                self.y_forces[beam_node_index_1]])
            beam.z_forces_global = np.array([
                self.z_forces[beam_node_index_0],
                self.z_forces[beam_node_index_1]])

            beam.x_moments_global = np.array([
                self.x_moments[beam_node_index_0],
                self.x_moments[beam_node_index_1]])
            beam.y_moments_global = np.array([
                self.y_moments[beam_node_index_0],
                self.y_moments[beam_node_index_1]])
            beam.z_moments_global = np.array([
                self.z_moments[beam_node_index_0],
                self.z_moments[beam_node_index_1]])

            beam.solve_shape_functions()


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

        # for shape functions
        eik_z = (self.material_properties["elastic modulus"] *
                 self.cross_section_properties["second moment of area z"] *
                 self.cross_section_properties["transverse shear deflection constant z"])

        eik_y = (self.material_properties["elastic modulus"] *
                 self.cross_section_properties["second moment of area y"] *
                 self.cross_section_properties["transverse shear deflection constant y"])

        ga = self.material_properties["shear modulus"] * self.cross_section_properties["area"]

        self.g_y = eik_y/ga
        self.g_z = eik_z/ga

    def cross_section_properties_init(self):
        if self.cross_section.lower() == "circle":
            return cross_section_circle(self.cross_section_parameters)
        elif self.cross_section.lower() == "annulus":
            return cross_section_annulus(self.cross_section_parameters)
        elif self.cross_section.lower() == "rectangle":
            return cross_section_rectangle(self.cross_section_parameters)
        elif self.cross_section.lower() == "hexagon":
            return cross_section_hexagon(self.cross_section_parameters)
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'

#%% TESTING
def test_cantilever_rectangle():
    beam_system = BeamSystem("Cantilever-end load")
    arm_beam = BeamType("rectangle", [20, 6], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([0,0,0]),"origin")
    beam_system.add_node(np.array([500,0,0]),"point_1")
    beam_system.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system.add_boundary_condition(0,"x","origin")
    beam_system.add_boundary_condition(0,"y","origin")
    beam_system.add_boundary_condition(0,"z","origin")
    beam_system.add_boundary_condition(0,"theta x","origin")
    beam_system.add_boundary_condition(0,"theta y","origin")
    beam_system.add_boundary_condition(0,"theta z","origin")
    beam_system.add_force(np.array([0, 0, -4]), "point_1")
    beam_system.add_boundary_condition(0,"x","point_1")
    beam_system.add_boundary_condition(0,"y","point_1")
    beam_system.add_boundary_condition(0,"theta x","point_1")
    beam_system.add_boundary_condition(0,"theta z","point_1")

    beam_system.solve_FEA()
    L = 500
    I = (6 * 20 ** 3) / 12
    E = 71700
    F = 4

    y = - (F * L ** 3) / (3*E*I)
    m = F * L
    dydx = (F*L**2)/(2*E*I) - (F*L**2)/(E*I)

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info("Started Cantilever—end load")
    print("calculated z displacement: " + str(y) + " | FEA z displacement: " + str(beam_system.z_displacements[1]))
    print("calculated y displacement: " + str(0) + " | FEA y displacement: " + str(beam_system.y_displacements[1]))
    print("calculated moment: " + str(-m) + " | FEA moment reaction: " + str(beam_system.y_moments[0]))
    print("calculated reaction: " + str(F) + " | FEA reaction z: " + str(beam_system.z_forces[0]))
    assert_almost_equal(y, beam_system.z_displacements[1], 3)
    assert_almost_equal(0, beam_system.y_displacements[1])
    assert_almost_equal(-dydx, beam_system.y_angles[1])
    assert_almost_equal(-m, beam_system.y_moments[0])
    assert_almost_equal(F, beam_system.z_forces[0])

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape z disp " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][0]))

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape y angle " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][1]))


def test_center_multidim():
    beam_system = BeamSystem("Multidimensional Simple Support Center Load")
    arm_beam = BeamType("annulus", [25, 20], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([-1000, -1000, -1000]), "point_0")
    beam_system.add_node(np.array([0, 0, 0]), "point_1")
    beam_system.add_node(np.array([1000, 1000, 1000]), "point_2")
    beam_system.add_beam(arm_beam, "point_0", "point_1", np.array([0, -1, 1]), "beam_1")
    beam_system.add_beam(arm_beam, "point_1", "point_2", np.array([0, -1, 1]), "beam_2")
    beam_system.add_boundary_condition(0, "x", "point_0")
    beam_system.add_boundary_condition(0, "y", "point_0")
    beam_system.add_boundary_condition(0, "z", "point_0")
    beam_system.add_boundary_condition(0, "x", "point_2")
    beam_system.add_boundary_condition(0, "y", "point_2")
    beam_system.add_boundary_condition(0, "z", "point_2")
    beam_system.add_force(np.array([0, 1, -1]), "point_1")
    beam_system.solve_FEA()

    L = math.sqrt(3 * 2000 ** 2)
    I = (math.pi/64) * (50**4 - 40**4)
    E = 71700
    F = math.sqrt(2)
    # kuu
    # Z1  Z2
    # Z2  Z3
    y = -(F * (L**3)) / (48 * E * I)
    r = F/2
    dydx = -(F * L ** 2) / (16 * E * I)
    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info('Started Multidimensional Simple Support Center Load Test')
    logging.debug("FEA x: " + str(beam_system.x_displacements[1]))
    logging.debug("FEA y: " + str(beam_system.y_displacements[1]))
    logging.debug("FEA z: " + str(beam_system.z_displacements[1]))
    logging.debug("FEA x theta: " + str(beam_system.x_angles[1]))
    logging.debug("FEA y theta: " + str(beam_system.y_angles[1]))
    logging.debug("FEA z theta: " + str(beam_system.z_angles[1]))
    logging.debug("FEA x reaction: " + str(beam_system.x_forces[0]))
    logging.debug("FEA y reaction: " + str(beam_system.y_forces[0]))
    logging.debug("FEA z reaction: " + str(beam_system.z_forces[0]))

    total_fea_disp = math.sqrt(beam_system.x_displacements[1] ** 2 + beam_system.y_displacements[1] ** 2 + beam_system.z_displacements[1] ** 2)
    total_fea_reaction = math.sqrt(beam_system.x_forces[0] ** 2 + beam_system.y_forces[0] ** 2 + beam_system.z_forces[0] ** 2)

    logging.debug("Calculated total displacement: " + str(-y) + " | Total displacement FEA: " + str(total_fea_disp))
    logging.debug("Calculated total reaction: " + str(r) + " | Total reaction FEA: " + str(total_fea_reaction))
    assert_almost_equal(-y, total_fea_disp, 3)
    assert_almost_equal(r, total_fea_reaction)
    logging.info('Finished')


def test_incline_boundary_conditions():
    beam_system = BeamSystem("Incline Boundary Conditions")
    beam_12 = BeamType("circle", [27.63953195], "Aluminum7075-T6", "beam_12")
    beam_3 = BeamType("circle", [32.869128059], "Aluminum7075-T6", "beam_3")

    beam_12.material_properties['elastic modulus'] = 210000
    beam_3.material_properties['elastic modulus'] = 210000
    beam_12.cross_section_properties['transverse shear deflection constant y'] = 0
    beam_12.cross_section_properties['transverse shear deflection constant z'] = 0
    beam_12.cross_section_properties['transverse shear deflection constant y'] = 0
    beam_12.cross_section_properties['transverse shear deflection constant z'] = 0
    print(beam_12.cross_section_properties)

    beam_system.add_node(np.array([0, 0, 0]), "point_1")
    beam_system.add_node(np.array([0, 1000, 0]), "point_2")
    beam_system.add_node(np.array([1000, 1000, 0]), "point_3")
    beam_system.add_beam(beam_12, "point_1", "point_2", np.array([0, 0, 1]), "beam_1")
    beam_system.add_beam(beam_12, "point_2", "point_3", np.array([0, 0, 1]), "beam_2")
    beam_system.add_beam(beam_3, "point_1", "point_3", np.array([0, 0, 1]), "beam_3")

    beam_system.add_boundary_condition(0, "x", "point_1")
    beam_system.add_boundary_condition(0, "y", "point_1")
    beam_system.add_boundary_condition(0, "z", "point_1")
    beam_system.add_boundary_condition(0, "theta x", "point_1")
    beam_system.add_boundary_condition(0, "theta y", "point_1")
    beam_system.add_boundary_condition(0, "theta z", "point_1")

    beam_system.add_boundary_condition(0, "y", "point_2")
    beam_system.add_boundary_condition(0, "z", "point_2")
    beam_system.add_boundary_condition(0, "theta y", "point_2")
    beam_system.add_boundary_condition(0, "theta z", "point_2")


    beam_system.add_boundary_condition(0, "y", "point_3", np.array([0.707, 0.707, 0]), np.array([-0.707, 0.707, 0]))
    beam_system.add_boundary_condition(0, "z", "point_3")
    beam_system.add_boundary_condition(0, "theta x", "point_3")
    beam_system.add_boundary_condition(0, "theta y", "point_3")
    beam_system.add_boundary_condition(0, "theta z", "point_3")
    beam_system.add_force(np.array([1000000, 0, 0]), "point_2")
    beam_system.solve_FEA()
    np.set_printoptions(linewidth=400)

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info('Started Incline Boundary Conditions')
    logging.debug("\nPRE-TRANSFORM 3:\n" + str(beam_system.global_stiffness_matrix))
    logging.debug("\nTRANSFORM 3:\n" + str(beam_system.local_bc_transform))
    logging.debug("\nPOST-TRANSFORM 3:\n" + str(beam_system.applied_stiffness_matrix))
    logging.debug("FEA 2y: " + str(beam_system.x_displacements[1]) + ' | Expected 11.91')
    # assert_almost_equal(beam_system.x_displacements[1], 11.91, 3)
    logging.debug("FEA 1x: " + str(beam_system.x_forces[0]) + ' | Expected -500000')
    logging.debug("FEA 1x: " + str(beam_system.x_forces[0]) + ' | Expected -500000')
    # assert_almost_equal(beam_system.x_forces[0], -500000, -4)
    logging.debug("FEA 1y: " + str(beam_system.y_forces[0]) + ' | Expected -500000')
    logging.debug("FEA 1y: " + str(beam_system.y_forces[0]) + ' | Expected -500000')
    # assert_almost_equal(beam_system.y_forces[0], -500000, -4)
    logging.debug("FEA 2y: " + str(beam_system.y_forces[1]) + ' | Expected 0')
    # assert_almost_equal(beam_system.y_forces[1], 0, -4)
    logging.info('Finished')


def test_reversed_cantilever_rectangle():
    L = 500
    I = (6 * 20 ** 3) / 12
    E = 71700
    F = 4

    y = - (F * L ** 3) / (3*E*I)
    m = F * L
    dydx = (F*L**2)/(2*E*I) - (F*L**2)/(E*I)

    beam_system = BeamSystem("Reverse Cantilever-end load")
    arm_beam = BeamType("rectangle", [20, 6], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([0,0,0]),"origin")
    beam_system.add_node(np.array([500,0,0]),"point_1")
    beam_system.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system.add_boundary_condition(0,"x","origin")
    beam_system.add_boundary_condition(0,"y","origin")
    beam_system.add_boundary_condition(0,"z","origin")
    beam_system.add_boundary_condition(0,"theta x","origin")
    beam_system.add_boundary_condition(0,"theta y","origin")
    beam_system.add_boundary_condition(0,"theta z","origin")
    beam_system.add_boundary_condition( -0.5816413736845836,"z", "point_1")


    beam_system.solve_FEA()

    print("\nFORCE MOMENT VECTOR: \n" + str(beam_system.force_moment_vector))

    print("KUU: " + str(beam_system.kuu))
    print("KPP: " + str(beam_system.kpp))

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info("Started Reverse Cantilever—end load")
    print("calculated z displacement: " + str(y) + " | FEA z displacement: " + str(beam_system.z_displacements[1]))
    print("calculated reaction: " + str(F) + " | FEA reaction z: " + str(beam_system.z_forces[1]))
    print("stress: ", beam_system.beams[0].stresses_bending_z)
    print("expected stress: ", m/I)

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape z disp " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][0]))

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape y angle " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][1]))

def test_cantilever_annulus():
    beam_system = BeamSystem("annulus Cantilever-end load")
    arm_beam = BeamType("hexagon", [20], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([0,0,0]),"origin")
    beam_system.add_node(np.array([500,0,0]),"point_1")
    beam_system.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system.add_boundary_condition(0,"x","origin")
    beam_system.add_boundary_condition(0,"y","origin")
    beam_system.add_boundary_condition(0,"z","origin")
    beam_system.add_boundary_condition(0,"theta x","origin")
    beam_system.add_boundary_condition(0,"theta y","origin")
    beam_system.add_boundary_condition(0,"theta z","origin")
    beam_system.add_force(np.array([0, 0, -4]), "point_1")

    beam_system.solve_FEA()

    beam_system_2 = BeamSystem("annulus Cantilever-end load")
    beam_system_2.add_node(np.array([0,0,0]),"origin")
    beam_system_2.add_node(np.array([500,0,0]),"point_1")
    beam_system_2.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system_2.add_boundary_condition(0,"x","origin")
    beam_system_2.add_boundary_condition(0,"y","origin")
    beam_system_2.add_boundary_condition(0,"z","origin")
    beam_system_2.add_boundary_condition(0,"theta x","origin")
    beam_system_2.add_boundary_condition(0,"theta y","origin")
    beam_system_2.add_boundary_condition(0,"theta z","origin")
    beam_system_2.add_boundary_condition(-0.00917048964187906,"z","point_1")

    beam_system_2.solve_FEA()
    L = 500
    I = (6 * 20 ** 3) / 12
    E = 71700
    F = 4

    y = - (F * L ** 3) / (3*E*I)
    m = F * L
    dydx = (F*L**2)/(2*E*I) - (F*L**2)/(E*I)

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info("Started annulus Cantilever—end load")
    print(" FEA z displacement: " + str(beam_system.z_displacements[1]))
    print(" FEA z force: " + str(beam_system_2.z_forces[1]))
    print("stress: ", beam_system.beams[0].stresses_bending_z)
    print("expected stress: ", m/beam_system.beams[0].beam_type.cross_section_properties['second moment of area y'])


    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape z disp " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][0]))

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape y angle " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][1]))

if __name__ == '__main__':
    test_cantilever_rectangle()
    test_center_multidim()
    test_incline_boundary_conditions()
    test_reverse_cantilever_rectangle()
    test_cantilever_annulus()

#%%
