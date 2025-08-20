import numpy as np

from Scripts.Finite_Elements import FEA_3D as FEA_3D


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
