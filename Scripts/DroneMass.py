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

        self.lumped_mass_frame = self.drone.blade_num * (self.arm_mass + self.strut_mass)

        self.arm_beam_centroids = np.zeros([self.drone.blade_num, 3])
        self.strut_beam_centroids = np.zeros([self.drone.blade_num, 3])
        self.centroid_calculation()
        print(self.lumped_mass_frame)
        print(self.strut_beam_centroids)

    def centroid_calculation(self):
        for i in range(self.drone.blade_num):
            theta = i * 2 * math.pi / self.drone.blade_num
            self.arm_beam_centroids[i][0] = math.cos(theta) * self.drone.nominal_rad / 2
            self.arm_beam_centroids[i][1] = math.sin(theta) * self.drone.nominal_rad / 2
            self.arm_beam_centroids[i][2] = 0

            self.strut_beam_centroids[i][0] = math.cos(theta + math.pi / self.drone.blade_num) * self.drone.strut_pos
            self.strut_beam_centroids[i][1] = math.sin(theta + math.pi / self.drone.blade_num) * self.drone.strut_pos
            self.strut_beam_centroids[i][2] = 0


if __name__ == '__main__':
    battery_1 = Battery(100, np.array([-5, 0, -1]), np.array([1, 1, 1]))
    battery_2 = Battery(100, np.array([5, 0, -1]), np.array([1, 1, 1]))

    arm_beam = BeamType("annulus", [25, 20], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [30, 20], "Aluminum7075-T6", "strut_beam")

    d1 = Drone("drone_test",500,200,8, arm_beam, strut_beam)
    
    d1_mass = DroneMass(d1, [battery_1, battery_2])

#%%
