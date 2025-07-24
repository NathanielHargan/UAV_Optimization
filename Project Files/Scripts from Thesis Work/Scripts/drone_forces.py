import numpy as np
from Scripts.material_properties import Gravity

air_density = 1.225e-12 # Mg/mm^3


class DroneForces:
    def __init__(self, mass, projected_surface_area, drag_coeff, profile):
        self.mass = mass
        self.projected_surface_area = projected_surface_area
        self.profile = profile
        self.drag_coeff = drag_coeff

    def time_solve(self, t):
        x_result, y_result, name = self.profile.time_solve(t)
        # x = x_result[0]
        # y = y_result[0]
        # v_x = x_result[1]
        v_y = y_result[1]
        # a_x = x_result[2]
        a_y = y_result[2]
        #j_x = x_result[3]
        # j_y = y_result[3]

        force_drag_y_mag = self.projected_surface_area * air_density * self.drag_coeff * (v_y ** 2) / 2
        force_drag_y_dir = 0
        if v_y > 0:
            force_drag_y_dir = -1
        else:
            force_drag_y_dir = 1

        force_drag_y = force_drag_y_mag * force_drag_y_dir

        # Gravity is negative
        # drag vector is the opposite direction as
        # a_y * self.mass = - Gravity * self.mass + force_drag_y + thrust_y

        thrust_y = (a_y * self.mass) + (Gravity * self.mass) - force_drag_y

        yield thrust_y
        yield force_drag_y

        return

    def time_solve_int(self, t):
        # C = 0
        v_y = y_result[1]
        a_y = y_result[2]

        force_drag_y_mag_int = self.projected_surface_area * air_density * self.drag_coeff * (v_y ** 2) / 2
        force_drag_y_dir = 0

        if v_y > 0:
            force_drag_y_dir = -1
        else:
            force_drag_y_dir = 1

        force_drag_y = force_drag_y_mag_int * force_drag_y_dir

        int_thrust_y = (a_y * self.mass) + (Gravity * self.mass) - force_drag_y
        return int_thrust_y


if __name__ == '__main__':
    from Scripts.drone_mass import DroneMass
    from Scripts.drone_battery import Battery
    from Scripts.drone_geometry import DroneGeometry
    from Scripts.beams import BeamType
    from Scripts.mission_profile import MissionProfileCubicSpline

    # Mass Analysis Setup
    battery_1 = Battery(4.292e-3, np.array([0, 86.75, -50]), np.array([260.5, 123.5, 63.5]))
    battery_2 = Battery(4.292e-3, np.array([0, -86.75, -50]), np.array([260.5, 123.5, 63.5]))
    arm_beam = BeamType("annulus", [29, 25], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [14, 12], "Aluminum7075-T6", "strut_beam")
    d1 = DroneGeometry("drone_test", 1000, 660, 8, arm_beam, strut_beam, [battery_1, battery_2])
    d1_mass = DroneMass(d1)

    # Mass Analysis results
    total_mass = d1_mass.total_mass
    surface_area = d1_mass.total_mass

    # Mission Profile generation
    m1 = MissionProfileCubicSpline("m1", 1, 1, 0)
    m1.add_segment("takeoff",10,50,5,5,1)
    m1.add_segment("climb",30,550,1000,30,1)
    m1.add_segment("cruise",100,1000,1100,30,1)
    m1.add_segment("descent",30,1100,500,0,0)
    m1.add_segment("loiter",30,1100,500,0,0)
    m1.add_segment("land",20,1200,0,0,0)

    # Drone Force Analysis
    d1_forces = DroneForces(total_mass, surface_area, 1, m1)
#%%
