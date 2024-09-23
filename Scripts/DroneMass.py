from Scripts.drone import Drone
from Scripts.beams import BeamType
import numpy as np

class Battery:
    def __init__(self, mass, x, y, z, l, h ,w):
        self.mass = mass
        self.position = np.array([x,y,z])
        self.dimensions = np.array([l,h,w])

        self.ixx = (1/12) * mass * (h**2 + w**2)
        self.iyy = (1/12) * mass * (l**2 + w**2)
        self.iyy = (1/12) * mass * (l**2 + h**2)

class DroneMass:
    def __init__(self, drone):
        self.drone = drone

        arm_linear_density = self.drone.arm_beam.material_properties["mass density"] * self.drone.arm_beam.cross_section_properties["area"]
        strut_linear_density = self.drone.strut_beam.material_properties["mass density"] * self.drone.strut_beam.cross_section_properties["area"]

        self.arm_mass = arm_linear_density * self.drone.nominal_rad
        self.strut_mass = strut_linear_density * self.drone.strut_length

        self.lumped_mass_frame = self.drone.blade_num * (self.arm_mass + self.strut_mass)

    def add_battery(self):



if __name__ == '__main__':
    arm_beam = BeamType("annulus", [25, 20], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [30, 20], "Aluminum7075-T6", "strut_beam")

    d1 = Drone("drone_test",500,200,8, arm_beam, strut_beam)
    
    d1_mass = DroneMass(d1)

#%%
