import numpy as np

class BeamMass:
    def __init__(self, beam_system):
        self.beams = beam_system.beams
        beam_num = len(self.beams)

        self.beam_mass = np.zeros(beam_num)
        self.beam_moment = np.zeros([beam_num, 3, 3])
        self.beam_centroids = np.zeros([beam_num, 3])
        self.total_mass = 0
        self.center_of_mass = np.zeros(3)
        self.beam_moment_of_inertia = np.zeros([beam_num, 3, 3])
        self.total_moment_of_inertia = np.zeros([3, 3])

        for i, beam in enumerate(self.beams):
            length = beam.length
            area = beam.beam_type.cross_section_properties["area"]
            density = beam.beam_type.material_properties["mass density"]
            mass = length * area * density

            self.beam_mass[i] = mass
            self.total_mass += mass

            transform = beam.transformation_matrix_pos
            self.beam_moment[i] =  transform @ beam.mass_moment_matrix @ np.transpose(transform)

            start_point = beam.start_node.location
            end_point = beam.end_node.location

            self.beam_centroids[i] = (start_point + end_point) / 2


        for i, beam_cent in enumerate(self.beam_centroids):
            self.center_of_mass += beam_cent * self.beam_mass[i] / self.total_mass

        for i, beam_cent in enumerate(self.beam_centroids):
            displacement = beam_cent - self.center_of_mass
            translation = self.beam_mass[i] * ((np.dot(displacement, displacement) * np.identity(3)) - np.outer(displacement, displacement)) #  Find source
            self.beam_moment_of_inertia[i] = self.beam_moment[i] + translation
            self.total_moment_of_inertia += self.beam_moment_of_inertia[i]

        eigenvalues, eigenvectors = np.linalg.eig(self.total_moment_of_inertia)
        self.total_moment_of_inertia_ev = eigenvalues
        self.total_moment_of_inertia_principal = np.sum(eigenvalues) * (1/3)



