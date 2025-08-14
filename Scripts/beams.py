#%%

#%%


#%%

#%%
import numpy as np
import math
import scipy
from Scripts.cross_section_properties import cross_section_circle
from Scripts.cross_section_properties import cross_section_annulus
from Scripts.cross_section_properties import cross_section_rectangle
from Scripts.cross_section_properties import cross_section_hexagon
from Scripts.cross_section_properties import cross_section_i_beam
from Scripts.material_properties import Materials
from Scripts.material_properties import NU_Daniel_failure
import Scripts.mass_properties as MassProperties
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
        self.transient_force = np.array([0, 0, 0, 0, 0, 0])

        self.dynamic_harmonic_force = np.zeros(6)
        self.dynamic_transient_force = np.empty(0)

        self.local_bc_transform = np.identity(6)
        self.result_displacement = np.array([0, 0, 0])
        self.result_angular_displacement = np.array([0, 0, 0])
        self.result_force = np.array([0, 0,  0])
        self.result_moment = np.array([0, 0, 0])

        self.result_displacement_amplitude = np.array([0, 0, 0])
        self.result_angular_displacement_amplitude = np.array([0, 0, 0])
        self.result_force_amplitude = np.array([0, 0, 0])
        self.result_moment_amplitude = np.array([0, 0, 0])


    def add_force(self, force):
        self.force = self.force + force

    def add_dynamic_harmonic_force(self, dynamic_force):
        self.dynamic_harmonic_force = self.dynamic_harmonic_force + dynamic_force

    def add_moment(self, moment):
        self.moment = self.moment + moment

    def add_dynamic_transient_force(self, dynamic_force_matrix):
        self.dynamic_transient_force = self.dynamic_transient_force + dynamic_force_matrix

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

        self.d_local_dynamic_harmonic = np.empty(12)
        self.r_local_dynamic_harmonic = np.empty(12)

        self.local_stiffness_matrix = np.empty((12, 12))
        self.local_mass_matrix = np.empty((12, 12))
        self.global_mass_matrix = np.empty((12, 12))
        self.local_bc_transformation = np.identity(12)
        self.transformation_matrix = np.empty((12, 12))
        self.transformation_matrix_node = np.empty((12, 12))
        self.transformation_matrix_node_inv = np.empty((12, 12))
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

        self.x_displacements_local_dynamic_harmonic = np.empty(2)
        self.y_displacements_local_dynamic_harmonic = np.empty(2)
        self.z_displacements_local_dynamic_harmonic = np.empty(2)
        self.x_angles_local_dynamic_harmonic = np.empty(2)
        self.y_angles_local_dynamic_harmonic = np.empty(2)
        self.z_angles_local_dynamic_harmonic = np.empty(2)
        self.x_forces_local_dynamic_harmonic = np.empty(2)
        self.y_forces_local_dynamic_harmonic = np.empty(2)
        self.z_forces_local_dynamic_harmonic = np.empty(2)
        self.x_moments_local_dynamic_harmonic = np.empty(2)
        self.y_moments_local_dynamic_harmonic = np.empty(2)
        self.z_moments_local_dynamic_harmonic = np.empty(2)

        self.x_displacements_local_ts = np.empty(2)
        self.y_displacements_local_ts = np.empty(2)
        self.z_displacements_local_ts = np.empty(2)
        self.x_angles_local_ts = np.empty(2)
        self.y_angles_local_ts = np.empty(2)
        self.z_angles_local_ts = np.empty(2)
        self.x_forces_local_ts = np.empty(2)
        self.y_forces_local_ts = np.empty(2)
        self.z_forces_local_ts = np.empty(2)
        self.x_moments_local_ts = np.empty(2)
        self.y_moments_local_ts = np.empty(2)
        self.z_moments_local_ts = np.empty(2)

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

        self.x_displacements_global_dynamic = np.empty(2)
        self.y_displacements_global_dynamic = np.empty(2)
        self.z_displacements_global_dynamic = np.empty(2)
        self.x_angles_global_dynamic = np.empty(2)
        self.y_angles_global_dynamic = np.empty(2)
        self.z_angles_global_dynamic = np.empty(2)
        self.x_forces_global_dynamic = np.empty(2)
        self.y_forces_global_dynamic = np.empty(2)
        self.z_forces_global_dynamic = np.empty(2)
        self.x_moments_global_dynamic = np.empty(2)
        self.y_moments_global_dynamic = np.empty(2)
        self.z_moments_global_dynamic = np.empty(2)

        self.x_displacements_global_ts = np.empty(2)
        self.y_displacements_global_ts = np.empty(2)
        self.z_displacements_global_ts = np.empty(2)
        self.x_angles_global_ts = np.empty(2)
        self.y_angles_global_ts = np.empty(2)
        self.z_angles_global_ts = np.empty(2)
        self.x_forces_global_ts = np.empty(2)
        self.y_forces_global_ts = np.empty(2)
        self.z_forces_global_ts = np.empty(2)
        self.x_moments_global_ts = np.empty(2)
        self.y_moments_global_ts = np.empty(2)
        self.z_moments_global_ts = np.empty(2)

        self.a_inv_y_shape_vector = np.empty(4)
        self.a_inv_z_shape_vector = np.empty(4)

        self.a_inv_y_shape_vector_ts = np.empty(4)
        self.a_inv_z_shape_vector_ts = np.empty(4)

        self.a_inv_y_shape_vector_dynamic_harmonic = np.empty(4)
        self.a_inv_z_shape_vector_dynamic_harmonic = np.empty(4)

        self.solve_stiffness_matrix()

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
        self.mass_moment_matrix = self.beam_type.mass_moment(self.length)
        self.local_mass_matrix = FEA_3D.mass_matrix_lumped_3d(
            self.beam_type.cross_section_properties["area"],
            self.length,
            self.beam_type.material_properties["mass density"],
            self.mass_moment_matrix
        )

        self.transformation_matrix = FEA_3D.transformation_matrix_element(self.direction, self.k_node)
        self.transformation_matrix_node = FEA_3D.transformation_matrix_node(self.direction, self.k_node)
        self.transformation_matrix_node_inv = np.linalg.inv(self.transformation_matrix_node)

        self.global_stiffness_matrix = FEA_3D.local_to_global_matrix(
            self.local_stiffness_matrix,
            self.transformation_matrix) #  np.transpose(t) @ k @ t

        self.global_mass_matrix = FEA_3D.local_to_global_matrix(self.local_mass_matrix, self.transformation_matrix)

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

    def solve_local_ts(self, d, r):
        self.d_local_ts = np.zeros(np.shape(d))
        self.r_local_ts = np.zeros(np.shape(r))
        for i in range(len(d)):
            self.d_local_ts[i] = self.transformation_matrix @ d[i]
            self.r_local_ts[i] = self.transformation_matrix @ r[i]

        self.x_displacements_local_ts = np.array([self.d_local_ts[:, 0], self.d_local_ts[:, 6]]).T
        self.y_displacements_local_ts = np.array([self.d_local_ts[:, 1], self.d_local_ts[:, 7]]).T
        self.z_displacements_local_ts = np.array([self.d_local_ts[:, 2], self.d_local_ts[:, 8]]).T

        self.x_angles_local_ts = np.array([self.d_local_ts[:, 3], self.d_local_ts[:, 9]]).T
        self.y_angles_local_ts = np.array([self.d_local_ts[:, 4], self.d_local_ts[:, 10]]).T
        self.z_angles_local_ts = np.array([self.d_local_ts[:, 5], self.d_local_ts[:, 11]]).T

        self.x_forces_local_ts = np.array([self.r_local_ts[:, 0], self.r_local_ts[:,6]]).T
        self.y_forces_local_ts = np.array([self.r_local_ts[:, 1], self.r_local_ts[:,7]]).T
        self.z_forces_local_ts = np.array([self.r_local_ts[:,2], self.r_local_ts[:,8]]).T

        self.x_moments_local_ts = np.array([self.r_local_ts[:,3], self.r_local_ts[:,9]]).T
        self.y_moments_local_ts = np.array([self.r_local_ts[:,4], self.r_local_ts[:,10]]).T
        self.z_moments_local_ts = np.array([self.r_local_ts[:,5], self.r_local_ts[:,11]]).T

    def solve_local_dynamic_harmonic(self, d, r):
        self.d_local_dynamic_harmonic = self.transformation_matrix @ d
        self.r_local_dynamic_harmonic = self.transformation_matrix @ r
        self.x_displacements_local_dynamic_harmonic = np.array([self.d_local_dynamic_harmonic[0], self.d_local_dynamic_harmonic[6]])
        self.y_displacements_local_dynamic_harmonic = np.array([self.d_local_dynamic_harmonic[1], self.d_local_dynamic_harmonic[7]])
        self.z_displacements_local_dynamic_harmonic = np.array([self.d_local_dynamic_harmonic[2], self.d_local_dynamic_harmonic[8]])

        self.x_angles_local_dynamic_harmonic = np.array([self.d_local_dynamic_harmonic[3], self.d_local_dynamic_harmonic[9]])
        self.y_angles_local_dynamic_harmonic = np.array([self.d_local_dynamic_harmonic[4], self.d_local_dynamic_harmonic[10]])
        self.z_angles_local_dynamic_harmonic = np.array([self.d_local_dynamic_harmonic[5], self.d_local_dynamic_harmonic[11]])

        self.x_forces_local_dynamic_harmonic = np.array([self.r_local_dynamic_harmonic[0], self.r_local_dynamic_harmonic[6]])
        self.y_forces_local_dynamic_harmonic = np.array([self.r_local_dynamic_harmonic[1], self.r_local_dynamic_harmonic[7]])
        self.z_forces_local_dynamic_harmonic = np.array([self.r_local_dynamic_harmonic[2], self.r_local_dynamic_harmonic[8]])
        self.x_moments_local_dynamic_harmonic = np.array([self.r_local_dynamic_harmonic[3], self.r_local_dynamic_harmonic[9]])
        self.y_moments_local_dynamic_harmonic = np.array([self.r_local_dynamic_harmonic[4], self.r_local_dynamic_harmonic[10]])
        self.z_moments_local_dynamic_harmonic = np.array([self.r_local_dynamic_harmonic[5], self.r_local_dynamic_harmonic[11]])

    def find_stresses(self):
        strn0, strss0, transv0, shearstrss0 = self.strain_stress_at_point_circle(0)
        strn1, strss1, transv1, shearstrss1 = self.strain_stress_at_point_circle(1)
        self.bending_strains_0 = strn0
        self.bending_strains_1 = strn1
        self.bending_stresses_0 = strss0
        self.bending_stresses_1 = strss1
        self.strain_transverse_shear_0 = transv0
        self.strain_transverse_shear_1 = transv1
        self.stress_transverse_shear_0 = shearstrss0
        self.stress_transverse_shear_1 = shearstrss1

    def find_stresses_dynamic(self):
        strn0, strss0, transv0, shearstrss0 = self.strain_stress_at_point_circle_dynamic(0)
        strn1, strss1, transv1, shearstrss1 = self.strain_stress_at_point_circle_dynamic(1)
        self.bending_strains_0_dynamic = strn0
        self.bending_strains_1_dynamic = strn1
        self.bending_stresses_0_dynamic = strss0
        self.bending_stresses_1_dynamic = strss1
        self.strain_transverse_shear_0_dynamic = transv0
        self.strain_transverse_shear_1_dynamic = transv1
        self.stress_transverse_shear_0_dynamic = shearstrss0
        self.stress_transverse_shear_1_dynamic = shearstrss1

    def find_stresses_ts(self):
        strn0, strss0, transv0, shearstrss0 = self.strain_stress_at_point_circle_ts(0)
        strn1, strss1, transv1, shearstrss1 = self.strain_stress_at_point_circle_ts(1)
        self.bending_strains_0_ts = strn0
        self.bending_strains_1_ts = strn1
        self.bending_stresses_0_ts = strss0
        self.bending_stresses_1_ts = strss1
        self.strain_transverse_shear_0_ts = transv0
        self.strain_transverse_shear_1_ts = transv1
        self.stress_transverse_shear_0_ts = shearstrss0
        self.stress_transverse_shear_1_ts = shearstrss1

    def strain_stress_at_point_circle(self, xi):
        y, z, dy, dz, ddy, ddz = self.return_shape_functions(xi)
        h = self.beam_type.cross_section_properties['circumscribed']
        strain_y, stress_y = self.stress_strain_shape_function(y, dy, h)
        strain_z, stress_z = self.stress_strain_shape_function(z, dz, h)

        strain_bending = np.linalg.norm([strain_y[0], strain_z[0]])
        stress_bending = np.linalg.norm([stress_y[0], stress_z[0]])
        strain_transverse_shear = strain_y[1]
        stress_transverse_shear = stress_y[1]

        return strain_bending, stress_bending, strain_transverse_shear, stress_transverse_shear

    def strain_stress_at_point_circle_dynamic(self, xi):
        y, z, dy, dz, ddy, ddz = self.return_shape_functions_dynamic(xi)
        h = self.beam_type.cross_section_properties['circumscribed']
        strain_y, stress_y = self.stress_strain_shape_function(y, dy, h)
        strain_z, stress_z = self.stress_strain_shape_function(z, dz, h)

        strain_bending = np.linalg.norm([strain_y[0], strain_z[0]])
        stress_bending = np.linalg.norm([stress_y[0], stress_z[0]])
        strain_transverse_shear = strain_y[1]
        stress_transverse_shear = stress_y[1]

        return strain_bending, stress_bending, strain_transverse_shear, stress_transverse_shear

    def strain_stress_at_point_circle_ts(self, xi):
        t_len = np.shape(self.x_displacements_global_ts)[0]
        strain_bending = np.zeros([t_len])
        stress_bending = np.zeros([t_len])
        strain_transverse_shear = np.zeros([t_len])
        stress_transverse_shear = np.zeros([t_len])
        h = self.beam_type.cross_section_properties['circumscribed']

        for ti in range(t_len):
            y, z, dy, dz, ddy, ddz = self.return_shape_functions_ts(xi, ti)
            strain_y, stress_y = self.stress_strain_shape_function(y, dy, h)
            strain_z, stress_z = self.stress_strain_shape_function(z, dz, h)

            strain_bending[ti] = np.linalg.norm([strain_y[0], strain_z[0]])
            stress_bending[ti] = np.linalg.norm([stress_y[0], stress_z[0]])
            strain_transverse_shear[ti] = strain_y[1]
            stress_transverse_shear[ti] = stress_y[1]

        return strain_bending, stress_bending, strain_transverse_shear, stress_transverse_shear

    def stress_strain_shape_function(self, n, dn, h):
        # Shape function[0] = uy(x)
        # Shape function[1] = thetaz(x)
        # Shape function[2] = gamma(x)
        # h is the height on the axis of the shape function (y or z local axis)
        strain_bending = -h * dn[1]
        stress_bending = strain_bending * self.beam_type.material_properties["elastic modulus"]
        strain_transverse_shear = n[2]
        stress_transverse_shear = self.beam_type.material_properties["shear modulus"] * n[2] / self.beam_type.cross_section_properties['transverse shear deflection constant y']
        strain_vector = np.array([strain_bending, strain_transverse_shear])
        stress_vector = np.array([stress_bending, stress_transverse_shear])
        return strain_vector, stress_vector

    def solve_shape_functions(self):
        a_inv_y = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_y)
        a_inv_z = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_z)

        self.a_inv_y_shape_vector = a_inv_y @ np.array([self.y_displacements_local[0],
                                                        self.z_angles_local[0],
                                                        self.y_displacements_local[1],
                                                        self.z_angles_local[1]])

        self.a_inv_z_shape_vector = a_inv_z @ np.array([self.z_displacements_local[0],
                                                        -self.y_angles_local[0],
                                                        self.z_displacements_local[1],
                                                        -self.y_angles_local[1]])

    def solve_shape_functions_ts(self):
        t_len = np.shape(self.x_displacements_global_ts)[0]

        a_inv_y = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_y)
        a_inv_z = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_z)

        self.a_inv_y_shape_vector_ts = np.zeros([t_len, 4])
        self.a_inv_z_shape_vector_ts = np.zeros([t_len, 4])

        for i in range(t_len):

            self.a_inv_y_shape_vector_ts[i] = a_inv_y @ np.array([self.y_displacements_local_ts[i,0],
                                                        self.z_angles_local_ts[i,0],
                                                        self.y_displacements_local_ts[i,1],
                                                        self.z_angles_local_ts[i,1]])

            self.a_inv_z_shape_vector_ts[i] = a_inv_z @ np.array([self.z_displacements_local_ts[i, 0],
                                                            -self.y_angles_local_ts[i, 0],
                                                            self.z_displacements_local_ts[i,1],
                                                            -self.y_angles_local_ts[i, 1]])

    def solve_shape_functions_dynamic_harmonic(self):
        a_inv_y = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_y)
        a_inv_z = FEA_3D.shape_function_timoshenko_a_inv(self.length, self.beam_type.g_z)

        self.a_inv_y_shape_vector_dynamic_harmonic = a_inv_y @ np.array([self.y_displacements_local_dynamic_harmonic[0],
                                                                         self.z_angles_local_dynamic_harmonic[0],
                                                                         self.y_displacements_local_dynamic_harmonic[1],
                                                                         self.z_angles_local_dynamic_harmonic[1]])

        self.a_inv_z_shape_vector_dynamic_harmonic = a_inv_z @ np.array([self.z_displacements_local_dynamic_harmonic[0],
                                                                         -self.y_angles_local_dynamic_harmonic[0],
                                                                         self.z_displacements_local_dynamic_harmonic[1],
                                                                         -self.y_angles_local_dynamic_harmonic[1]])

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

        # derivative with respect to x (both y and z are the same because g is unchanging with respect to x)

        dX_y_z = np.array([
            [0, 1, 2 * x, 2 * x**2],
            [0, 0, 2, (6 * x)],
            [0, 0, 0, 0]])

        ddX_y_z = np.array([
            [0, 0, 2, 4 * x],
            [0, 0, 0, 6],
            [0, 0, 0, 0]])

        y = X_y @ self.a_inv_y_shape_vector
        z = X_z @ self.a_inv_z_shape_vector
        dy = dX_y_z @ self.a_inv_y_shape_vector
        dz = dX_y_z @ self.a_inv_z_shape_vector
        ddy = ddX_y_z @ self.a_inv_y_shape_vector
        ddz = ddX_y_z @ self.a_inv_z_shape_vector
        return y, z, dy, dz, ddy, ddz

    def return_shape_functions_dynamic(self, xi):
        x = xi * self.length
        X_y = np.array([
            [1, x, x**2, x**3],
            [0, 1, 2*x, (3 * x ** 2) - 6 * self.beam_type.g_y],
            [0, 0, 0, -6 * self.beam_type.g_y]])

        X_z = np.array([
            [1, x, x**2, x**3],
            [0, 1, 2*x, (3 * x ** 2) - 6 * self.beam_type.g_z],
            [0, 0, 0, -6 * self.beam_type.g_z]])

        # derivative with respect to x (both y and z are the same because g is unchanging with respect to x)

        dX_y_z = np.array([
            [0, 1, 2 * x, 2 * x**2],
            [0, 0, 2, (6 * x)],
            [0, 0, 0, 0]])

        ddX_y_z = np.array([
            [0, 0, 2, 4 * x],
            [0, 0, 0, 6],
            [0, 0, 0, 0]])

        y = X_y @ self.a_inv_y_shape_vector_dynamic_harmonic
        z = X_z @ self.a_inv_z_shape_vector_dynamic_harmonic
        dy = dX_y_z @ self.a_inv_y_shape_vector_dynamic_harmonic
        dz = dX_y_z @ self.a_inv_z_shape_vector_dynamic_harmonic
        ddy = ddX_y_z @ self.a_inv_y_shape_vector_dynamic_harmonic
        ddz = ddX_y_z @ self.a_inv_z_shape_vector_dynamic_harmonic
        return y, z, dy, dz, ddy, ddz

    def return_shape_functions_ts(self, xi, ti):
        x = xi * self.length
        X_y = np.array([
            [1, x, x**2, x**3],
            [0, 1, 2*x, (3 * x ** 2) - 6 * self.beam_type.g_y],
            [0, 0, 0, -6 * self.beam_type.g_y]])

        X_z = np.array([
            [1, x, x**2, x**3],
            [0, 1, 2*x, (3 * x ** 2) - 6 * self.beam_type.g_z],
            [0, 0, 0, -6 * self.beam_type.g_z]])

        # derivative with respect to x (both y and z are the same because g is unchanging with respect to x)

        dX_y_z = np.array([
            [0, 1, 2 * x, 2 * x**2],
            [0, 0, 2, (6 * x)],
            [0, 0, 0, 0]])

        ddX_y_z = np.array([
            [0, 0, 2, 4 * x],
            [0, 0, 0, 6],
            [0, 0, 0, 0]])

        y = X_y @ self.a_inv_y_shape_vector_ts[ti]
        z = X_z @ self.a_inv_z_shape_vector_ts[ti]
        dy = dX_y_z @ self.a_inv_y_shape_vector_ts[ti]
        dz = dX_y_z @ self.a_inv_z_shape_vector_ts[ti]
        ddy = ddX_y_z @ self.a_inv_y_shape_vector_ts[ti]
        ddz = ddX_y_z @ self.a_inv_z_shape_vector_ts[ti]

        return y, z, dy, dz, ddy, ddz

    def return_deformation_global(self, xi):
        n_y, n_z, dn_y, dn_z, ddn_y, ddn_z = self.return_shape_functions(xi)
        n_x = self.x_displacements_local[0] + xi * (self.x_displacements_local[1] - self.x_displacements_local[0]) # linear
        n_x_theta = self.x_angles_local[0] + xi * (self.x_angles_local[1] - self.x_angles_local[0]) # linear
        deformation_local = np.array([[n_x], [n_y[0]], [n_z[0]], [n_x_theta], [n_y[1]], [n_z[1]]])
        deformation_global = self.transformation_matrix_node_inv @ deformation_local
        return deformation_global

    def failure_criterion(self):
        max_shear = max(self.stress_transverse_shear_0, self.stress_transverse_shear_1)
        max_stress = max(self.bending_stresses_0, self.bending_stresses_1)
        return NU_Daniel_failure(max_stress, max_shear, self.beam_type.material_properties)

    def calc_failure_criterion_swt(self):
        fatigue_str_coeff = self.beam_type.material_properties["fatigue strength coefficient"]
        fatigue_str_exp = self.beam_type.material_properties["fatigue strength exponent"]
        self.stress_max_0 = self.bending_stresses_0 + self.bending_stresses_0_dynamic
        self.stress_max_1 = self.bending_stresses_1 + self.bending_stresses_1_dynamic
        self.stress_equivalent_alternating_0 = np.sqrt(self.stress_max_0 * self.bending_stresses_0_dynamic)
        self.stress_equivalent_alternating_1 = np.sqrt(self.stress_max_1 * self.bending_stresses_1_dynamic)
        self.cycle_num_0 = (1/2) * (self.stress_equivalent_alternating_0/fatigue_str_coeff) ** (1/fatigue_str_exp)
        self.cycle_num_1 = (1/2) * (self.stress_equivalent_alternating_1/fatigue_str_coeff) ** (1/fatigue_str_exp)

    def failure_criterion_swt(self):

        cycle_limit = self.beam_type.material_properties["fatigue strength cycles"]

        return np.array([self.cycle_num_0 / cycle_limit, self.cycle_num_1 / cycle_limit])



class BeamSystem:
    def __init__(self, name):
        self.system_name = name  # A string with the name of the system
        self.beam_node_indexes = np.empty((0, 2))  # List of the node indexes that beams go between

        self.boundary_conditions_node_indexes = np.empty((0, 1))
        self.boundary_conditions_type = np.empty((0, 1))
        self.boundary_conditions = np.empty((0, 1))

        self.initial_conditions_node_indexes_position = np.empty((0, 1))
        self.initial_conditions_type_position = np.empty((0, 1))
        self.initial_conditions_position = np.empty((0, 1))

        self.initial_conditions_node_indexes_velocity = np.empty((0, 1))
        self.initial_conditions_type_velocity = np.empty((0, 1))
        self.initial_conditions_velocity = np.empty((0, 1))

        self.initial_conditions_node_indexes_acceleration = np.empty((0, 1))
        self.initial_conditions_type_acceleration = np.empty((0, 1))
        self.initial_conditions_acceleration = np.empty((0, 1))

        self.nodes = []  # list of node objects
        self.beams = []  # List of beam objects

        # list of node names for calling
        self.node_names = []  # list of node names
        self.beam_names = []  # list of beam names

        # placeholder

        self.global_length = 0

        self.global_stiffness_matrix = np.empty(0)
        self.global_mass_matrix = np.empty(0)
        self.bc_transform = np.empty(0)
        self.applied_stiffness_matrix = np.empty(0)
        self.applied_mass_matrix = np.empty(0)

        self.ru = np.empty(0)
        self.dp = np.empty(0)

        self.p_index = np.empty(0)
        self.u_index = np.empty(0)

        self.displacement_angle_vector = np.empty(0)
        self.force_moment_vector = np.empty(0)

        self.up_index = np.empty(0)

        self.kuu = np.empty(0)
        self.kup = np.empty(0)
        self.kpu = np.empty(0)
        self.kpp = np.empty(0)

        self.muu = np.empty(0)
        self.mup = np.empty(0)
        self.mpu = np.empty(0)
        self.mpp = np.empty(0)

        self.du = np.empty(0)
        self.dup = np.empty(0)
        self.ru = np.empty(0)
        self.rup = np.empty(0)

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

        self.mag_displacements = np.empty(0)

        self.frequency = np.empty(0)

        self.x_displacements_amplitude = np.empty(0)
        self.y_displacements_amplitude = np.empty(0)
        self.z_displacements_amplitude = np.empty(0)

        self.x_angles_amplitude = np.empty(0)
        self.y_angles_amplitude = np.empty(0)
        self.z_angles_amplitude = np.empty(0)

        self.x_forces_amplitude = np.empty(0)
        self.y_forces_amplitude = np.empty(0)
        self.z_forces_amplitude = np.empty(0)

        self.x_moments_ts= np.empty(0)
        self.y_moments_ts = np.empty(0)
        self.z_moments_ts = np.empty(0)

        self.x_displacements_ts = np.empty(0)
        self.y_displacements_ts = np.empty(0)
        self.z_displacements_ts = np.empty(0)

        self.x_angles_ts = np.empty(0)
        self.y_angles_ts = np.empty(0)
        self.z_angles_ts = np.empty(0)

        self.x_forces_ts = np.empty(0)
        self.y_forces_ts = np.empty(0)
        self.z_forces_ts = np.empty(0)

        self.x_moments_ts = np.empty(0)
        self.y_moments_ts = np.empty(0)
        self.z_moments_ts = np.empty(0)

        self.a_int_constants = np.empty(0)

        self.kdtuu = np.empty(0)
        self.kdtup = np.empty(0)
        self.kdtpp = np.empty(0)
        self.kdtpu = np.empty(0)

        self.timesteps = np.empty(0)
        self.timestep_num = 0

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
            return [bc_ref]
        elif type(bc_ref) is str:
            return {
                "x": [0],
                "y": [1],
                "z": [2],
                "theta x": [3],
                "theta y": [4],
                "theta z": [5],
                "x symm": [0, 4, 5],
                "y symm": [1, 3, 5],
                "z symm": [2, 3, 4],
                "pinned": [0, 1, 2],
                "en castre": [0, 1, 2, 3, 4, 5],
            }[bc_ref]  # life hack
        else:
            print("invalid boundary condition type")
            return "error"

    def add_boundary_condition(self, boundary_val, bc_type_ref, node_ref, direction=None, k_node_dir=None):
        node_index = self.select_node(node_ref)
        bc_types = self.select_boundary_condition_type(bc_type_ref)
        for bc_type in bc_types:
            self.boundary_conditions = np.concatenate((self.boundary_conditions, [[boundary_val]]))
            self.boundary_conditions_type = np.concatenate((self.boundary_conditions_type, [[bc_type]]))
            self.boundary_conditions_node_indexes = np.concatenate((self.boundary_conditions_node_indexes, [[node_index]]))

            if direction is not None:
                node = self.nodes[node_index]
                node.create_local_transform(direction, k_node_dir)

    def add_initial_condition_position(self, initial_val, initial_type_ref, node_ref, direction=None, k_node_dir=None):
        node_index = self.select_node(node_ref)
        initial_types = self.select_boundary_condition_type(initial_type_ref)

        for initial_type in initial_types:
            self.initial_conditions_position = np.concatenate((self.initial_conditions_position, [[initial_val]]))
            self.initial_conditions_type_position = np.concatenate((self.initial_conditions_type_position, [[initial_type]]))
            self.initial_conditions_node_indexes_position = np.concatenate((self.initial_conditions_node_indexes_position, [[node_index]]))

            if direction is not None:
                node = self.nodes[node_index]
                node.create_local_transform(direction, k_node_dir)

    def add_initial_condition_velocity(self, initial_val, initial_type_ref, node_ref, direction=None, k_node_dir=None):
        node_index = self.select_node(node_ref)
        initial_types = self.select_boundary_condition_type(initial_type_ref)

        for initial_type in initial_types:
            self.initial_conditions_velocity = np.concatenate((self.initial_conditions_velocity, [[initial_val]]))
            self.initial_conditions_type_velocity = np.concatenate((self.initial_conditions_type_velocity, [[initial_type]]))
            self.initial_conditions_node_indexes_velocity = np.concatenate((self.initial_conditions_node_indexes_velocity, [[node_index]]))

            if direction is not None:
                node = self.nodes[node_index]
                node.create_local_transform(direction, k_node_dir)

    def add_initial_condition_acceleration(self, initial_val, initial_type_ref, node_ref, direction=None, k_node_dir=None):
        node_index = self.select_node(node_ref)
        initial_types = self.select_boundary_condition_type(initial_type_ref)

        for initial_type in initial_types:
            self.initial_conditions_acceleration = np.concatenate((self.initial_conditions_acceleration, [[initial_val]]))
            self.initial_conditions_type_acceleration = np.concatenate((self.initial_conditions_type_acceleration, [[initial_type]]))
            self.initial_conditions_node_indexes_acceleration = np.concatenate((self.initial_conditions_node_indexes_acceleration, [[node_index]]))

            if direction is not None:
                node = self.nodes[node_index]
                node.create_local_transform(direction, k_node_dir)


    def add_force(self, direction, node_ref):
        node = self.nodes[self.select_node(node_ref)]
        node.add_force(direction)

    def init_dynamic_forces(self, frequency):
        self.frequency = frequency

    def add_dynamic_harmonic_force(self, direction, node_ref):
        node = self.nodes[self.select_node(node_ref)]
        node.add_dynamic_harmonic_force(direction)

    def add_dynamic_transient_force(self, direction_func, node_ref):
        direction_matrix = np.zeros([self.timestep_num, 6])
        for i, t in enumerate(self.timesteps):
            direction_matrix[i] = direction_func(t)
        node = self.nodes[self.select_node(node_ref)]
        node.add_dynamic_transient_force(direction_matrix)

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
        self.global_length = np.shape(self.nodes)[0] * 6

        dp = self.boundary_conditions[:, 0]
        p_index = np.intc(self.boundary_conditions_node_indexes * 6 + self.boundary_conditions_type)

        u_index_predel = np.intc(np.arange(self.global_length))
        ru_predel = np.zeros(self.global_length)  # fu before deleting bcs

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

        self.up_index = np.concatenate((u_index, p_index), axis=None)
        global_element_matrices = np.empty((len(self.beam_node_indexes), 12, 12))
        global_element_mass_matrices = np.empty((len(self.beam_node_indexes), 12, 12))


        for i in range(len(self.beam_node_indexes)):
            global_element_matrices[i] = self.beams[i].global_stiffness_matrix
            global_element_mass_matrices[i] = self.beams[i].global_mass_matrix

        # assembles the global bc node transforms
        self.bc_transform = np.zeros([len(self.nodes) * 6, len(self.nodes) * 6])
        for i, node in enumerate(self.nodes):
            self.bc_transform[i*6:i*6 + 6,i*6:i*6 + 6] = node.local_bc_transform

        self.global_stiffness_matrix = self.bc_transform @ FEA_3D.assemble_stiffness_3d(self.beam_node_indexes, global_element_matrices, self.global_length) @ np.transpose(self.bc_transform)
        self.global_mass_matrix = self.bc_transform @ FEA_3D.assemble_stiffness_3d(self.beam_node_indexes, global_element_mass_matrices, self.global_length) @ np.transpose(self.bc_transform)

        self.applied_stiffness_matrix = FEA_3D.rearrange_stiffness_matrix(self.global_stiffness_matrix, self.up_index)
        self.applied_mass_matrix = FEA_3D.rearrange_stiffness_matrix(self.global_mass_matrix, self.up_index)

        kuu, kup, kpu, kpp = FEA_3D.partition_stiffness_matrix(self.applied_stiffness_matrix, len(u_index))
        muu, mup, mpu, mpp = FEA_3D.partition_stiffness_matrix(self.applied_mass_matrix, len(u_index))

        self.kuu = kuu
        self.kup = kup
        self.kpu = kpu
        self.kpp = kpp

        self.muu = muu
        self.mup = mup
        self.mpu = mpu
        self.mpp = mpp

        self.displacement_angle_vector_amplitude = np.empty([len(self.nodes)*6])
        self.force_moment_vector_amplitude = np.empty([len(self.nodes)*6])

    def solve_static(self):
        self.du = scipy.linalg.solve(self.kuu, self.ru - (self.kup @ self.dp))
        self.rp = (self.kpu @ self.du) + (self.kpp @ self.dp)

        dup = np.concatenate((self.du, self.dp), axis=None)
        rup = np.concatenate((self.ru, self.rp), axis=None)

        self.rup = rup
        self.dup = dup

        d = np.empty(self.global_length) # d'
        r = np.empty(self.global_length) # r'

        for i in range(self.global_length):
            d[self.up_index[i]] = dup[i]
            r[self.up_index[i]] = rup[i]

        self.displacement_angle_vector = np.transpose(self.bc_transform) @ d
        self.force_moment_vector = np.transpose(self.bc_transform) @ r

        self.x_displacements = self.displacement_angle_vector[0::6]
        self.y_displacements = self.displacement_angle_vector[1::6]
        self.z_displacements = self.displacement_angle_vector[2::6]
        self.mag_displacements = np.linalg.norm(np.vstack((self.x_displacements, self.y_displacements, self.z_displacements)),axis = 0)

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
            local_d = np.concatenate((self.displacement_angle_vector[beam_index_0:beam_index_0 + 6], self.displacement_angle_vector[beam_index_1:beam_index_1+6]))
            local_r = np.concatenate((self.force_moment_vector[beam_index_0:beam_index_0 + 6], self.force_moment_vector[beam_index_1:beam_index_1+6]))
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
            beam.find_stresses()

    def solve_failure(self):
        max_failure_crit = 0
        for beam in self.beams:
            fail_crit = beam.failure_criterion()
            max_failure_crit = max([fail_crit[0], fail_crit[1], max_failure_crit])
        self.max_failure_crit = max_failure_crit

    def solve_failure_swt(self):
        max_failure_crit_swt = 0
        for beam in self.beams:
            fail_crit_swt = beam.failure_criterion_swt()
            max_failure_crit_swt = max([fail_crit_swt[0], fail_crit_swt[1], max_failure_crit_swt])
        self.max_failure_crit_swt = max_failure_crit_swt

    def solve_natural_frequencies(self):
        eigenvals = scipy.linalg.eigvals(self.kuu, self.muu)
        self.natural_frequencies_eigenvals = eigenvals
        self.natural_frequencies = np.sort(np.real(np.sqrt(eigenvals)))

    def solve_dynamic_harmonic(self, alpha, beta):
        self.displacement_angle_vector_amplitude = np.empty([self.global_length])
        self.force_moment_vector_amplitude = np.empty([self.global_length])

        ru_predel = np.zeros(self.global_length)  # fu before deleting bcs

        for i in range(len(self.nodes)):  # for loop + counter
            r_index = int(i * 6)
            ru_predel[r_index:r_index+6] += self.nodes[i].dynamic_harmonic_force

        ru = np.delete(ru_predel, self.p_index, axis=0) # deleting boundary conditions

        ang_vel = self.frequency

        self.kduu = self.kuu - (ang_vel ** 2) * self.muu - ang_vel * (self.muu * alpha + self.kuu * beta)
        self.kdup = self.kup - (ang_vel ** 2) * self.mup - ang_vel * (self.mup * alpha + self.kup * beta)
        self.kdpp = self.kpp - (ang_vel ** 2) * self.mpp - ang_vel * (self.mpp * alpha + self.kpp * beta)
        self.kdpu = self.kpu - (ang_vel ** 2) * self.mpu - ang_vel * (self.mpu * alpha + self.kpu * beta)

        du = scipy.linalg.solve(self.kduu, ru - (self.kdup @ self.dp))
        rp = (self.kdpu @ du) + (self.kdpp @ self.dp)

        dup = np.concatenate((du, self.dp), axis=None)
        rup = np.concatenate((self.ru, rp), axis=None)

        d = np.empty(self.global_length)  # d'
        r = np.empty(self.global_length)  # r'

        for i in range(self.global_length):
            d[self.up_index[i]] = dup[i]
            r[self.up_index[i]] = rup[i]

        displacement_angle_vector = np.transpose(self.bc_transform) @ d
        force_moment_vector = np.transpose(self.bc_transform) @ r

        self.displacement_angle_vector_amplitude = displacement_angle_vector
        self.force_moment_vector_amplitude = force_moment_vector

        # Replace with Vstack
        self.x_displacements_amplitude = displacement_angle_vector[0::6]
        self.y_displacements_amplitude = displacement_angle_vector[1::6]
        self.z_displacements_amplitude = displacement_angle_vector[2::6]

        self.x_angles_amplitude = displacement_angle_vector[3::6]
        self.y_angles_amplitude = displacement_angle_vector[4::6]
        self.z_angles_amplitude = displacement_angle_vector[5::6]

        self.x_forces_amplitude = force_moment_vector[0::6]
        self.y_forces_amplitude = force_moment_vector[1::6]
        self.z_forces_amplitude = force_moment_vector[2::6]

        self.x_moments_amplitude = force_moment_vector[3::6]
        self.y_moments_amplitude = force_moment_vector[4::6]
        self.z_moments_amplitude = force_moment_vector[5::6]

        for beam_num, beam in enumerate(self.beams):
            beam_index_0 = int(self.beam_node_indexes[beam_num][0] * 6)
            beam_index_1 = int(self.beam_node_indexes[beam_num][1] * 6)
            local_d = np.concatenate((self.displacement_angle_vector_amplitude[beam_index_0:beam_index_0 + 6], self.displacement_angle_vector_amplitude[beam_index_1:beam_index_1+6]))
            local_r = np.concatenate((self.force_moment_vector_amplitude[beam_index_0:beam_index_0 + 6], self.force_moment_vector_amplitude[beam_index_1:beam_index_1+6]))
            beam.solve_local_dynamic_harmonic(local_d, local_r)

            beam_node_index_0 = int(self.beam_node_indexes[beam_num][0])
            beam_node_index_1 = int(self.beam_node_indexes[beam_num][1])

            beam.x_displacements_global_dynamic_harmonic = np.array([
                self.x_displacements_amplitude[beam_node_index_0],
                self.x_displacements_amplitude[beam_node_index_1]])
            beam.y_displacements_global_dynamic_harmonic = np.array([
                self.y_displacements_amplitude[beam_node_index_0],
                self.y_displacements_amplitude[beam_node_index_1]])
            beam.z_displacements_global_dynamic_harmonic = np.array([
                self.z_displacements_amplitude[beam_node_index_0],
                self.z_displacements_amplitude[beam_node_index_1]])

            beam.x_angles_global_dynamic_harmonic = np.array([
                self.x_angles_amplitude[beam_node_index_0],
                self.x_angles_amplitude[beam_node_index_1]])
            beam.y_angles_global_dynamic_harmonic = np.array([
                self.y_angles_amplitude[beam_node_index_0],
                self.y_angles_amplitude[beam_node_index_1]])
            beam.z_angles_global_dynamic_harmonic = np.array([
                self.z_angles_amplitude[beam_node_index_0],
                self.z_angles_amplitude[beam_node_index_1]])

            beam.x_forces_global_dynamic_harmonic = np.array([
                self.x_forces_amplitude[beam_node_index_0],
                self.x_forces_amplitude[beam_node_index_1]])
            beam.y_forces_global_dynamic_harmonic = np.array([
                self.y_forces_amplitude[beam_node_index_0],
                self.y_forces_amplitude[beam_node_index_1]])
            beam.z_forces_global_dynamic_harmonic = np.array([
                self.z_forces_amplitude[beam_node_index_0],
                self.z_forces_amplitude[beam_node_index_1]])

            beam.x_moments_global_dynamic_harmonic = np.array([
                self.x_moments_amplitude[beam_node_index_0],
                self.x_moments_amplitude[beam_node_index_1]])
            beam.y_moments_global_dynamic_harmonic = np.array([
                self.y_moments_amplitude[beam_node_index_0],
                self.y_moments_amplitude[beam_node_index_1]])
            beam.z_moments_global_dynamic_harmonic = np.array([
                self.z_moments_amplitude[beam_node_index_0],
                self.z_moments_amplitude[beam_node_index_1]])

            beam.solve_shape_functions_dynamic_harmonic()
            beam.find_stresses_dynamic()
            beam.calc_failure_criterion_swt()

    def initialize_dynamic_transient(self, dt, ts_num, alpha, delta):
        # dt = timestep interval
        # ts_num = number of timesteps

        # alpha and delta integration constants
        if delta < 0.5:
            print("delta should be greater than 0.5")

        if alpha < 0.25 * (delta + 0.5) ** 2:
            print("alpha should be greater than 0.25 * (delta + 0.5)^2")

        a0 = 1 / (alpha * dt**2)
        a1 = delta / (alpha * dt)
        a2 = 1 / (alpha * dt)
        a3 = 1 / (2 * alpha) - 1
        a4 = delta / alpha - 1
        a5 = (dt/2) * ((delta / alpha) - 2)
        a6 = dt * (1 - delta)
        a7 = delta * dt

        self.a_int_constants = np.array([a0, a1, a2, a3, a4, a5, a6, a7])

        self.timestep_num = ts_num
        self.timesteps = np.linspace(0, ts_num * dt, ts_num)

        for node in self.nodes:
            node.dynamic_transient_force = np.zeros([ts_num, 6])

    def solve_dynamic_transient(self):
        # unchanged
        # dp = self.dp

        a0, a1, a2, a3, a4, a5, a6, a7 = self.a_int_constants

        u_len = len(self.ru)
        p_len = len(self.dp)

        cuu = np.zeros([u_len, u_len]) # empty

        self.kdtuu = self.kuu + a0 * self.muu
        self.kdtup = self.kup + a0 * self.mup
        self.kdtpp = self.kpp + a0 * self.mpp
        self.kdtpu = self.kpu + a0 * self.mpu

        du_initial = np.zeros(u_len)
        du_d_initial = np.zeros(u_len)
        du_dd_initial = np.zeros(u_len)

        du = du_initial
        du_d = du_d_initial
        du_dd = du_dd_initial

        # Array of displacement vectors
        du_ts = np.zeros([self.timestep_num, u_len])
        du_d_ts = np.zeros([self.timestep_num, u_len])
        du_dd_ts = np.zeros([self.timestep_num, u_len])

        rp_ts = np.zeros([self.timestep_num, p_len])
        ru_ts = np.zeros([self.timestep_num, u_len])

        du_ts[0] = du_initial
        du_d_ts[0] = du_d_initial
        du_dd_ts[0] = du_dd_initial

        # Solve for U.
        for i, t in enumerate(self.timesteps[0:]):
            ru_predel = np.zeros(self.global_length)  # fu before deleting bcs
            for j in range(len(self.nodes)):  # for loop + counter
                r_index = int(j * 6)
                ru_predel[r_index:r_index+6] += self.nodes[j].dynamic_transient_force[i]

            ru = np.delete(ru_predel, self.p_index, axis=0)  # deleting boundary conditions

            ru_eff = ru + self.muu @ (a0 * du + a2 * du_d + a3 * du_dd)  # effective loads

            # Do include the kdtup @ dp term? It's not included in the matlab
            du_new = scipy.linalg.solve(self.kdtuu, ru_eff)

            # solving for RP?

            du_dd_new = a0 * (du_new - du) - a2 * du_d - a3 * du_dd

            du_d_new = du_d + a6 * du_dd + a7 * du_dd_new

            # update

            du = du_new
            du_d = du_d_new
            du_dd = du_dd_new

            ru_ts[i] = ru
            du_ts[i] = du
            du_d_ts[i] = du_d
            du_dd_ts[i] = du_dd

        # reaction forces
        for i, t in enumerate(self.timesteps):
            rp_ts[i] = self.mpu @ du_dd_ts[i] + self.kpu @ du_ts[i] + self.kpp @ self.dp

        dp_ts = np.tile(self.dp, (self.timestep_num, 1))
        dp_d_ts = np.zeros([self.timestep_num, p_len])
        dp_dd_ts = np.zeros([self.timestep_num, p_len])



        dup_ts = np.concatenate((du_ts, dp_ts), axis=1)
        dup_d_ts = np.concatenate((du_d_ts, dp_d_ts), axis=1)
        dup_dd_ts = np.concatenate((du_dd_ts, dp_dd_ts), axis=1)

        rup_ts = np.concatenate((ru_ts, rp_ts), axis=1)

        d_ts = np.empty([self.timestep_num, self.global_length])  # d
        d_d_ts = np.empty([self.timestep_num, self.global_length])  # d'
        d_dd_ts = np.empty([self.timestep_num, self.global_length])  # d''

        r_ts = np.empty([self.timestep_num, self.global_length])  # r'

        for i, t in enumerate(self.timesteps):
            for j in range(self.global_length):
                d_ts[i, self.up_index[j]] = dup_ts[i, j]
                d_d_ts[i, self.up_index[j]] = dup_d_ts[i, j]
                d_dd_ts[i, self.up_index[j]] = dup_dd_ts[i, j]
                r_ts[i, self.up_index[j]] = rup_ts[i,j]

        # transformation matrix
        for i, t in enumerate(self.timesteps):
            d_ts[i] = np.transpose(self.bc_transform) @ d_ts[i]
            d_d_ts[i] = np.transpose(self.bc_transform) @ d_d_ts[i]
            d_dd_ts[i] = np.transpose(self.bc_transform) @ d_dd_ts[i]
            r_ts[i] = np.transpose(self.bc_transform) @ r_ts[i]

        self.d_ts = d_ts
        self.d_d_ts = d_d_ts
        self.d_dd_ts = d_dd_ts

        self.r_ts = r_ts

        self.x_displacements_ts = d_ts[:,0::6]
        self.y_displacements_ts = d_ts[:,1::6]
        self.z_displacements_ts = d_ts[:,2::6]

        self.x_angles_ts = d_ts[:,3::6]
        self.y_angles_ts = d_ts[:,4::6]
        self.z_angles_ts = d_ts[:,5::6]

        self.x_forces_ts = r_ts[:,0::6]
        self.y_forces_ts = r_ts[:,1::6]
        self.z_forces_ts = r_ts[:,2::6]

        self.x_moments_ts = r_ts[:,3::6]
        self.y_moments_ts = r_ts[:,4::6]
        self.z_moments_ts = r_ts[:,5::6]

        for beam_num, beam in enumerate(self.beams):
            beam_index_0 = int(self.beam_node_indexes[beam_num][0] * 6)
            beam_index_1 = int(self.beam_node_indexes[beam_num][1] * 6)

            local_d = np.concatenate((self.d_ts[:, beam_index_0:beam_index_0 + 6], self.d_ts[:, beam_index_1:beam_index_1 + 6]), axis=1)
            local_r = np.concatenate((self.r_ts[:, beam_index_0:beam_index_0 + 6], self.r_ts[:, beam_index_1:beam_index_1 + 6]), axis=1)

            beam.solve_local_ts(local_d, local_r)

            beam_node_index_0 = int(self.beam_node_indexes[beam_num][0])
            beam_node_index_1 = int(self.beam_node_indexes[beam_num][1])

            beam.x_displacements_global_ts = np.array([
                self.x_displacements_ts[:, beam_node_index_0],
                self.x_displacements_ts[:, beam_node_index_1]]).T
            beam.y_displacements_global_ts = np.array([
                self.y_displacements_ts[:, beam_node_index_0],
                self.y_displacements_ts[:, beam_node_index_1]]).T
            beam.z_displacements_global_ts = np.array([
                self.z_displacements_ts[:, beam_node_index_0],
                self.z_displacements_ts[:, beam_node_index_1]]).T

            beam.x_angles_global_ts = np.array([
                self.x_angles_ts[:, beam_node_index_0],
                self.x_angles_ts[:, beam_node_index_1]]).T
            beam.y_angles_global_ts = np.array([
                self.y_angles_ts[:, beam_node_index_0],
                self.y_angles_ts[:, beam_node_index_1]]).T
            beam.z_angles_global_ts = np.array([
                self.z_angles_ts[:, beam_node_index_0],
                self.z_angles_ts[:, beam_node_index_1]]).T

            beam.x_forces_global_ts = np.array([
                self.x_forces_ts[:, beam_node_index_0],
                self.x_forces_ts[:, beam_node_index_1]]).T
            beam.y_forces_global_ts = np.array([
                self.y_forces_ts[:, beam_node_index_0],
                self.y_forces_ts[:, beam_node_index_1]]).T
            beam.z_forces_global_ts = np.array([
                self.z_forces_ts[:, beam_node_index_0],
                self.z_forces_ts[:, beam_node_index_1]]).T

            beam.x_moments_global_ts = np.array([
                self.x_moments_ts[:, beam_node_index_0],
                self.x_moments_ts[:, beam_node_index_1]]).T
            beam.y_moments_global_ts = np.array([
                self.y_moments_ts[:, beam_node_index_0],
                self.y_moments_ts[:, beam_node_index_1]]).T
            beam.z_moments_global_ts = np.array([
                self.z_moments_ts[:, beam_node_index_0],
                self.z_moments_ts[:, beam_node_index_1]]).T

            beam.solve_shape_functions_ts()
            beam.find_stresses_ts()



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



    def mass_moment(self, L):
        if self.cross_section.lower() == "circle":
            return MassProperties.mass_moment_circle(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "annulus":
            return MassProperties.mass_moment_annulus(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "rectangle":
            return MassProperties.mass_moment_rectangle(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "hexagon":
            return MassProperties.mass_moment_hexagon(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "i_beam":
            return MassProperties.mass_moment_i_beam(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'


    def cross_section_properties_init(self):
        if self.cross_section.lower() == "circle":
            return cross_section_circle(self.cross_section_parameters)
        elif self.cross_section.lower() == "annulus":
            return cross_section_annulus(self.cross_section_parameters)
        elif self.cross_section.lower() == "rectangle":
            return cross_section_rectangle(self.cross_section_parameters)
        elif self.cross_section.lower() == "hexagon":
            return cross_section_hexagon(self.cross_section_parameters)
        elif self.cross_section.lower() == "i_beam":
            return cross_section_i_beam(self.cross_section_parameters)
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'


def cantilever_rectangle():
    # Relevant
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
    beam_system.solve_static()
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


def center_multidim():
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


def incline_boundary_conditions():
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


def reversed_cantilever_rectangle():
    # Relevant
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


def cantilever_annulus():
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


def beam_30_deg():
    beam_system = BeamSystem("30 Deg Axial Tension")
    arm_beam = BeamType("circle", [20], "Aluminum7075-T6", "arm_beam")

    # 2 length 200 elements
    beam_system.add_node(np.array([0, 0, 0]),"origin")
    beam_system.add_node(np.array([math.sqrt(3)*100, 100, 0]),"point_1")
    beam_system.add_node(np.array([500, 200, 0]),"point_2")

    beam_system.add_beam(arm_beam, "origin", "point_1", np.array([0, 0, 1]), "axial_tension_beam")
    beam_system.add_beam(arm_beam, "point_1", "point_2", np.array([0, 0, 1]), "pulling_beam")

    beam_system.add_boundary_condition(0, "x", "origin")
    beam_system.add_boundary_condition(0, "y", "origin")
    beam_system.add_boundary_condition(0, "z", "origin")
    beam_system.add_boundary_condition(0, "theta x", "origin")
    beam_system.add_boundary_condition(0, "theta y", "origin")
    beam_system.add_boundary_condition(0, "theta z", "origin")

    beam_system.add_boundary_condition(0, "y", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "z", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "theta x", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "theta y", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "theta z", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))

    beam_system.add_force(np.array([0, 1, 0]), "point_2")

    beam_system.solve_FEA()
    beam_system.solve_static()

    N1_disp = beam_system.nodes[beam_system.select_node("point_1")].result_displacement
    print("Displacement 1: ", N1_disp)
    print("Displacement 2: ", beam_system.nodes[beam_system.select_node("point_2")].result_displacement)

    print("Displacement 1 theta: ", math.degrees(math.atan2(N1_disp[0], N1_disp[1])))


def timoshenko_simply_supp_point():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    midpoint_index = int((node_count-1)/2)
    L = 1000
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x+0.0001, 0, 0]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 1, 0]))

    # Create BCS
    beam_system.add_boundary_condition(0, "y", 0)
    beam_system.add_boundary_condition(0, "y", node_count-1)


    print("~Point Load Slender~")
    beam_system.add_force(np.array([0, -Q, 0]), midpoint_index)

    beam_system.solve_FEA()
    xi_midpoint = midpoint_index
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[midpoint_index].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -0.178522mm
    # ABAQUS 23 midpoint: -0.178871mm
    # Abaqus max stress: 1.920

    for beam in beam_system.beams:
        print("Max Stress 0", max(beam.bending_stresses_0))
        print("Max Stress 1", max(beam.bending_stresses_1))


def timoshenko_simply_supp():
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")

    node_count = 21
    L = 1000
    Q = 0.0222411 # Load
    E = arm_beam.material_properties['elastic modulus']
    I = arm_beam.cross_section_properties['second moment of area z']
    EI = E*I
    K = arm_beam.cross_section_properties['transverse shear deflection constant z']
    G = arm_beam.material_properties["shear modulus"]
    A = arm_beam.cross_section_properties['area']

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x+0.0001, 0, 0]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 1, 0]))

    # Create BCS
    beam_system.add_boundary_condition(0, "x", 0)
    beam_system.add_boundary_condition(0, "y", 0)
    beam_system.add_boundary_condition(0, "x", node_count-1)
    beam_system.add_boundary_condition(0, "y", node_count-1)

    for i in range(node_count):
        beam_system.add_boundary_condition(0, "z symm", i)
    '''
    beam_system.add_force(np.array([0, -Q, 0]), 10)
    beam_system.solve_FEA()

    print("Midpoint Disp:", beam_system.nodes[10].result_displacement)
    print("Midpoint Rotation:", beam_system.nodes[10].result_displacement)
    # Abaqus Results
    #
    '''

    # Add distributed load
    q = L*Q/(node_count)  # load on each node
    for i in range(1, node_count-1):  # ingnoring phantom loads
        beam_system.add_force(np.array([0, -q, 0]), i)
        print("")

    beam_system.add_moment(np.array([0, 0, (Q * (L/(node_count-1)) ** 2)/12]), 0)
    beam_system.add_moment(np.array([0, 0, -(Q * (L/(node_count-1)) ** 2)/12]), node_count-1)

    beam_system.solve_FEA()
    xi_midpoint = int(node_count/2)
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[xi_midpoint].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    x = beam_system.nodes[xi_midpoint].location[0]
    print("x:", x)
    w_bending = (Q/EI) * ((x**4 / 24) - ((L**2 * x**2)/16) + (5 * L ** 4)/384)
    w_shear = (Q/(2 * K * G * A)) * (((L ** 2) / 4) - (x ** 2))
    w = w_bending + w_shear
    print("Expected Midpoint Disp from bending: ", w_bending)
    print("Expected Midpoint Disp from shear: ", w_shear)
    print("Expected Midpoint Disp from total: ", w)


    x = beam_system.nodes[xi_75].location[0]
    print("75:",x)
    w_bending = (Q/EI) * ((x**4 / 24) - ((L**2 * x**2)/16) + (5 * L ** 4)/384)
    w_shear = (Q/(2 * K * G * A)) * (((L ** 2) / 4) - (x ** 2))
    w = w_bending + w_shear
    print("Expected 75 Disp from bending: ", w_bending)
    print("Expected 75 Disp from shear: ", w_shear)
    print("Expected 75 Disp from total: ", w)


def timoshenko_simply_supp_point_shear():
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    L = 50
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x, 0, 0]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 1, 0]))

    # Create BCS
    beam_system.add_boundary_condition(0, "x", 0)
    beam_system.add_boundary_condition(0, "x", node_count-1)
    beam_system.add_boundary_condition(0, "y", 0)
    beam_system.add_boundary_condition(0, "y", node_count-1)

    for i in range(node_count):
        beam_system.add_boundary_condition(0, "z symm", i)

    print("~Point Load Shear Test~")
    beam_system.add_force(np.array([0, -Q, 0]), 10)

    beam_system.solve_FEA()

    xi_midpoint = int(node_count/2)
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[10].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -3.997e-05
    # ABAQUS 23 midpoint: -2.232e-05


def timoshenko_stress_circle():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("circle", [25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    midpoint_index = int((node_count-1)/2)
    L = 1000
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x+0.0001, 0, 0.01]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    # Create BCS
    beam_system.add_boundary_condition(0, "pinned", 0)
    beam_system.add_boundary_condition(0, "pinned", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", 0)


    print("~Point Load Slender~")
    beam_system.add_force(np.array([0, -Q, -Q]), midpoint_index)

    beam_system.solve_FEA()
    xi_midpoint = midpoint_index
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[midpoint_index].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -0.178522mm
    # ABAQUS 23 midpoint: -0.178871mm
    # Abaqus max stress: 1.920
    for beam in beam_system.beams:
        print("Max Stress 0", np.round(beam.bending_stresses_0,5))
        print("Max Stress 1", np.round(beam.bending_stresses_1,5))

    print(beam_system.y_angles)
    # print(beam_system.du)

    # print(beam_system.kuu)
    print(beam_system.z_angles)
    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.du)


def timoshenko_stress():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25,25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    midpoint_index = int((node_count-1)/2)
    L = 1000
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x + 0.0001, 0.0001, 0.01]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    # Create BCS
    beam_system.add_boundary_condition(0, "pinned", 0)
    beam_system.add_boundary_condition(0, "pinned", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", 0)


    print("~Point Load Slender~")
    beam_system.add_force(np.array([0, -Q, -Q]), midpoint_index)

    beam_system.solve_FEA()
    xi_midpoint = midpoint_index
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[midpoint_index].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -0.178522mm
    # ABAQUS 23 midpoint: -0.178871mm
    # Abaqus max stress: 1.920
    for beam in beam_system.beams:
        print("Max Stress 0", np.round(beam.bending_stresses_0,5))
        print("Max Stress 1", np.round(beam.bending_stresses_1,5))

    # print(beam_system.y_angles)
    # print(beam_system.du)

    # print(beam_system.kuu)
    # print(beam_system.z_angles)
    print(np.shape(beam_system.global_stiffness_matrix))
    print(np.shape(beam_system.global_mass_matrix))

    beam_system.solve_natural_frequencies()
    print("nf:", beam_system.natural_frequencies)
    print("ev:", beam_system.natural_frequencies_eigenvals)
    print("ev shape:", np.shape(beam_system.natural_frequencies_eigenvals))
    print("k shape:", np.shape(beam_system.global_stiffness_matrix))

    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.du)


def timoshenko_freq():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25,25], "Aluminum7075-T6", "arm_beam")
    node_count = 50
    midpoint_index = int((node_count-1)/2)
    L = 1000

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x + 0.0001, 0.0001, 0.01]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    # Create BCS
    beam_system.add_boundary_condition(0, "pinned", 0)
    beam_system.add_boundary_condition(0, "pinned", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", 0)
    beam_system.solve_FEA()

    np.set_printoptions(linewidth=4000)
    # print(np.shape(beam_system.kuu))
    # print(np.shape(beam_system.muu))

    # print("kuu,", np.round(beam_system.kuu,7))
    print("muu,", np.round(beam_system.muu,7))
    # print("m,", np.round(beam_system.global_mass_matrix,7))
    # print("m_e,", np.round(beam_system.beams[0].global_mass_matrix,7))

    beam_system.solve_natural_frequencies()
    # print("lamda:", beam_system.natural_frequencies_eigenvals)
    # print("omega:", beam_system.natural_frequencies)
    print("Hz:", np.round(beam_system.natural_frequencies/(2*math.pi)))
    # print("ev shape:", np.shape(beam_system.natural_frequencies_eigenvals))
    # print("k shape:", np.shape(beam_system.global_stiffness_matrix))

    # m = 8.0605e-4 * np.array([[2,0],[0,1]])
    # k = 2.6507e6 * np.array([[2,-1],[-1,1]])
    # print("m:",m)
    # print("k:",k)
    # eigenvals = scipy.linalg.eigvals(k, m)
    # print("e,", eigenvals)
    # print("w,", np.sqrt(eigenvals))
    # print("hz,", np.sqrt(eigenvals)/(2*math.pi))


def transient_displacement():
    timesteps = 100
    node_count = 3
    L = 100

    beam_system = BeamSystem("canilever beam")
    arm_beam = BeamType("annulus", [14, 12], "Carbon Fiber", "arm_beam")

    for i in range(node_count):
        x = L * (i/(node_count-1))
        beam_system.add_node(np.array([x + 0.0001, 0.0001, 0.01]))

    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    beam_system.add_force(np.array([0,0,0.5]), node_count-1)
    beam_system.add_boundary_condition(0, "en castre", 0)

    def trans_force(t):
        return np.array([0, 0, (t * (t < 0.5) + (1 - t) * (t > 0.5)) * (t < 1), 0, 0, 0])

    beam_system.solve_FEA()
    beam_system.solve_static()
    beam_system.initialize_dynamic_transient(0.02, timesteps, 0.25,0.5)
    beam_system.add_dynamic_transient_force(trans_force, node_count-1)
    beam_system.solve_dynamic_transient()

    print("F(t) = ")
    print(beam_system.nodes[node_count-1].dynamic_transient_force)

    print("U_tip(t) = ")
    print(beam_system.d_ts[:, 2 + 6*(node_count - 1)])

    print("U_d_tip(t) = ")
    print(beam_system.d_d_ts[:, 2 + 6*(node_count - 1)])

    print("Stress(t) = ")
    print(beam_system.beams[0].bending_stresses_0_ts)

    print("U_static = ")
    print(beam_system.z_displacements[node_count-1])

    print("U_tip(0.5s) = ")
    print(beam_system.d_ts[int(0.5/0.02), 2 + 6*(node_count - 1)])



if __name__ == '__main__':
    # beam_30_deg()
    # timoshenko_simply_supp_point()
    # timoshenko_simply_supp()
    # timoshenko_simply_supp()
    # timoshenko_simply_supp_point_shear()
    # timoshenko_stress()
    # timoshenko_stress_circle()
    # timoshenko_freq()
    transient_displacement()
