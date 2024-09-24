from Scripts.drone import Drone
from Scripts.beams import BeamType
import numpy as np
import math

class Battery:
    def __init__(self, mass, position, dimensions):
        self.mass = mass
        self.position = position
        self.dimensions = dimensions

        l = dimensions[0]
        h = dimensions[1]
        w = dimensions[2]

        self.ixx = (1/12) * mass * (h**2 + w**2)
        self.iyy = (1/12) * mass * (l**2 + w**2)
        self.iyy = (1/12) * mass * (l**2 + h**2)

class DroneMass:
    def __init__(self, drone, batteries):
        self.drone = drone
        self.batteries = batteries

        arm_linear_density = self.drone.arm_beam.material_properties["mass density"] * self.drone.arm_beam.cross_section_properties["area"]
        strut_linear_density = self.drone.strut_beam.material_properties["mass density"] * self.drone.strut_beam.cross_section_properties["area"]

        self.arm_mass = arm_linear_density * self.drone.nominal_rad
        self.strut_mass = strut_linear_density * self.drone.strut_length

        self.arm_ixx = self.drone.arm_beam

        self.lumped_mass_frame = self.drone.blade_num * (self.arm_mass + self.strut_mass)

        self.arm_beam_centroids = np.zeros([self.drone.blade_num, 3])
        self.strut_beam_centroids = np.zeros([self.drone.blade_num, 3])
        self.battery_centroids = np.zeros([len(self.batteries), 3])
        self.centroid_calculation()

        self.total_mass = self.lumped_mass_frame

        self.total_mass_calculation() # add battery mass

        self.centroids = np.empty([])
        self.masses = np.empty([])
        self.center_of_mass = np.zeros(3)
        self.center_of_mass_calcuation()

    def centroid_calculation(self):
        for i in range(self.drone.blade_num):
            theta = i * 2 * math.pi / self.drone.blade_num
            self.arm_beam_centroids[i][0] = math.cos(theta) * self.drone.nominal_rad / 2
            self.arm_beam_centroids[i][1] = math.sin(theta) * self.drone.nominal_rad / 2
            self.arm_beam_centroids[i][2] = 0

            self.strut_beam_centroids[i][0] = math.cos(theta + math.pi / self.drone.blade_num) * self.drone.strut_pos
            self.strut_beam_centroids[i][1] = math.sin(theta + math.pi / self.drone.blade_num) * self.drone.strut_pos
            self.strut_beam_centroids[i][2] = 0

        for i, battery in enumerate(self.batteries):
            self.battery_centroids[i] = battery.position


    def total_mass_calculation(self):
        for battery in self.batteries:
            self.total_mass += battery.mass

    def center_of_mass_calcuation(self):
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
        print(self.center_of_mass)




if __name__ == '__main__':
    battery_1 = Battery(100, np.array([-5, 0, -1]), np.array([1, 1, 1]))
    battery_2 = Battery(100, np.array([5, 0, -1]), np.array([1, 1, 1]))

    arm_beam = BeamType("annulus", [25, 20], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [30, 20], "Aluminum7075-T6", "strut_beam")

    d1 = Drone("drone_test",500,200,8, arm_beam, strut_beam)
    
    d1_mass = DroneMass(d1, [battery_1, battery_2])

#%%
