import numpy as np
import scipy

from Scripts.Finite_Elements import FEA_3D as FEA_3D
from Scripts.Finite_Elements.Beam_FEA_System.beam_element import BeamElement
from Scripts.Finite_Elements.Beam_FEA_System.node import Node


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
        self.mag_displacements_amplitude = np.empty(0)
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

        new_beam = BeamElement(beam_type, self.nodes[start_node_index], self.nodes[end_node_index], k_node, name)
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

        self.mag_displacements_amplitude = np.linalg.norm(np.vstack((self.x_displacements_amplitude, self.y_displacements_amplitude, self.z_displacements_amplitude)),axis = 0)

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
            # beam.calc_failure_criterion_swt()

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
