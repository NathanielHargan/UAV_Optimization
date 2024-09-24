from Scripts.drone import Drone
from Scripts.beams import BeamType
import numpy as np
import math


class Battery:
    def __init__(self, mass, position, dimensions):
        self.mass = mass
        self.position = position # position of the center of the battery
        self.dimensions = dimensions

        l = dimensions[0]
        h = dimensions[1]
        w = dimensions[2]

        self.local_moment_of_inertia = np.array([[(1/12) * mass * (h**2 + w**2), 0, 0],
                                                [0, (1/12) * mass * (l**2 + w**2), 0],
                                                [0, 0, (1/12) * mass * (l**2 + h**2)]])

        self.transformation = np.array([[1, 0 ,0],
                                        [0, 1, 0],
                                        [0, 0, 1]])


class DroneMass:
    def __init__(self, drone, batteries):
        self.drone = drone
        self.batteries = batteries

        arm_linear_density = self.drone.arm_beam.material_properties["mass density"] * self.drone.arm_beam.cross_section_properties["area"]
        strut_linear_density = self.drone.strut_beam.material_properties["mass density"] * self.drone.strut_beam.cross_section_properties["area"]

        self.arm_mass = arm_linear_density * self.drone.nominal_rad
        self.strut_mass = strut_linear_density * self.drone.strut_length
        self.lumped_mass_frame = self.drone.blade_num * (self.arm_mass + self.strut_mass)

        self.total_mass = self.lumped_mass_frame
        self.total_mass_calculation()  # add battery masses

        # moment of inertia of a slender rod
        self.arm_local_moment_of_inertia = np.array([[(1/12) * self.arm_mass * self.drone.nominal_rad ** 2, 0, 0],
                                                     [0, 0, 0],
                                                     [0, 0, (1/12) * self.arm_mass * self.drone.nominal_rad ** 2]])

        self.strut_local_moment_of_inertia = np.array([[(1/12) * self.strut_mass * self.drone.nominal_rad ** 2, 0, 0],
                                                     [0, 0, 0],
                                                     [0, 0, (1/12) * self.strut_mass * self.drone.nominal_rad ** 2]])

        self.arm_beam_centroids = np.zeros([self.drone.blade_num, 3])
        self.strut_beam_centroids = np.zeros([self.drone.blade_num, 3])
        self.battery_centroids = np.zeros([len(self.batteries), 3])
        self.centroid_calculation()

        self.centroids = np.empty([])
        self.masses = np.empty([])
        self.center_of_mass = np.zeros(3)
        self.center_of_mass_calculation()

        self.battery_transformation = np.zeros([len(self.batteries), 3, 3])
        self.arm_transformation = np.zeros([self.drone.blade_num, 3, 3])
        self.strut_transformation = np.zeros([self.drone.blade_num, 3, 3])
        self.calculate_transformation_matrices()

        self.battery_global_moments_of_inertia = np.zeros([len(self.batteries), 3, 3])
        self.arm_global_moments_of_inertia = np.zeros([self.drone.blade_num, 3, 3])
        self.strut_global_moments_of_inertia = np.zeros([self.drone.blade_num, 3, 3])
        self.total_moment_of_inertia = np.zeros([3,3])
        self.calculate_global_moments()

    def centroid_calculation(self):
        for i in range(self.drone.blade_num):
            theta = i * 2 * math.pi / self.drone.blade_num
            self.arm_beam_centroids[i][0] = math.cos(theta) * self.drone.nominal_rad / 2
            self.arm_beam_centroids[i][1] = math.sin(theta) * self.drone.nominal_rad / 2
            self.arm_beam_centroids[i][2] = 0

            self.strut_beam_centroids[i][0] = (math.cos(theta) + math.cos(theta + 2 * math.pi / self.drone.blade_num)) * self.drone.strut_pos / 2
            self.strut_beam_centroids[i][1] = (math.sin(theta) + math.sin(theta + 2 * math.pi / self.drone.blade_num)) * self.drone.strut_pos / 2
            self.strut_beam_centroids[i][2] = 0

        for i, battery in enumerate(self.batteries):
            self.battery_centroids[i] = battery.position

    def total_mass_calculation(self):
        for battery in self.batteries:
            self.total_mass += battery.mass

    def transformation_matrix(self, coord_dir, k_node_dir):
        coord_unit = coord_dir / np.linalg.norm(coord_dir)
        k_node_unit = k_node_dir / np.linalg.norm(k_node_dir)
        ortho_dir = np.cross(coord_dir, k_node_dir)
        ortho_unit = ortho_dir / np.linalg.norm(ortho_dir)

        l1 = coord_unit[0]
        l2 = k_node_unit[0]
        l3 = ortho_unit[0]
        m1 = coord_unit[1]
        m2 = k_node_unit[1]
        m3 = ortho_unit[1]
        n1 = coord_unit[2]
        n2 = k_node_unit[2]
        n3 = ortho_unit[2]

        t = np.array([[l1, m1, n1],
                      [l2, m2, n2],
                      [l3, m3, n3]])
        return t

    def calculate_transformation_matrices(self):
        for i, battery in enumerate(self.batteries):
            self.battery_transformation[i] = battery.transformation

        for i in range(self.drone.blade_num):
            theta = i * 2 * math.pi / self.drone.blade_num
            self.arm_transformation[i] = self.transformation_matrix(
                np.array([math.cos(theta), math.sin(theta), 0]),
                np.array([0, 0, 1]))

            theta_minus_1 = (i-1) * 2 * math.pi / self.drone.blade_num
            self.strut_transformation[i] = self.transformation_matrix(
                np.array([math.cos(theta) - math.cos(theta_minus_1), math.sin(theta) - math.sin(theta_minus_1), 0]),
                np.array([0, 0, 1]))

    def center_of_mass_calculation(self):
        self.centroids = np.transpose(np.vstack([self.arm_beam_centroids, self.strut_beam_centroids, self.battery_centroids]))
        arm_beam_masses = np.zeros([self.drone.blade_num])
        strut_beam_masses = np.zeros([self.drone.blade_num])
        battery_masses = np.zeros([len(self.batteries)])

        for i in range(self.drone.blade_num):
            arm_beam_masses[i] = self.arm_mass
            strut_beam_masses[i] = self.strut_mass

        for i, battery in enumerate(self.batteries):
            battery_masses[i] = battery.mass

        self.masses = np.transpose(np.hstack([arm_beam_masses, strut_beam_masses, battery_masses]))
        self.center_of_mass = (self.centroids @ self.masses) / self.total_mass


    def calculate_global_moments(self):
        for i, battery in enumerate(self.batteries):
            rotated_moment = self.battery_transformation[i] @ self.batteries[i].local_moment_of_inertia @ np.transpose(self.battery_transformation[i])
            displacement = self.batteries[i].position - self.center_of_mass
            translation = self.batteries[i].mass * (np.dot(displacement, displacement) * np.identity(3) - np.outer(displacement, displacement))
            # https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
            self.battery_global_moments_of_inertia[i] = rotated_moment + translation
            self.total_moment_of_inertia += self.battery_global_moments_of_inertia[i]

        for i in range(self.drone.blade_num):
            rotated_moment = self.arm_transformation[i] @ self.arm_local_moment_of_inertia @ np.transpose(self.arm_transformation[i])
            displacement = self.arm_beam_centroids[i] - self.center_of_mass
            translation = self.arm_mass * (np.dot(displacement, displacement) * np.identity(3) - np.outer(displacement, displacement))
            self.arm_global_moments_of_inertia[i] = rotated_moment + translation
            self.total_moment_of_inertia += self.arm_global_moments_of_inertia[i]

        for i in range(self.drone.blade_num):
            rotated_moment = self.strut_transformation[i] @ self.strut_local_moment_of_inertia @ np.transpose(self.strut_transformation[i])
            displacement = self.strut_beam_centroids[i] - self.center_of_mass
            translation = self.strut_mass * (np.dot(displacement, displacement) * np.identity(3) - np.outer(displacement, displacement))
            self.strut_global_moments_of_inertia[i] = rotated_moment + translation
            self.total_moment_of_inertia += self.strut_global_moments_of_inertia[i]


if __name__ == '__main__':
    # 40.05 N weight
    # Mass in N-s^2/mm (Mg)
    battery_1 = Battery(4.292e-3, np.array([0, 86.75, -25]), np.array([260.5, 123.5, 63.5]))
    battery_2 = Battery(4.292e-3, np.array([0, -86.75, -25]), np.array([260.5, 123.5, 63.5]))

    # al: 0.098 lb/in^3
    # 26.6 kN/m^3 =* 10^-9 *
    arm_beam = BeamType("annulus", [29, 25], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [14, 12], "Aluminum7075-T6", "strut_beam")

    d1 = Drone("drone_test", 500, 200, 8, arm_beam, strut_beam)
    
    d1_mass = DroneMass(d1, [battery_1, battery_2])
    print("total mass: ", d1_mass.total_mass)
    print("moment of inertia tensor: ", d1_mass.total_moment_of_inertia)
    print("center of mass: ", d1_mass.center_of_mass)

    # TEST IN SOLIDWORKS

#%%
