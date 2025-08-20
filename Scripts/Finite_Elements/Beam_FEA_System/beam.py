import numpy as np

from Scripts.Finite_Elements import FEA_3D as FEA_3D
from Scripts.Finite_Elements.material_properties import NU_Daniel_failure


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
