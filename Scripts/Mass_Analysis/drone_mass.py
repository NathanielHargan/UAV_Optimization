import numpy as np
import math



class DroneMass:
    def __init__(self, drone):
        self.drone = drone
        self.batteries = drone.batteries

        arm_linear_density = self.drone.arm_beam.material_properties["mass density"] * self.drone.arm_beam.cross_section_properties["area"]
        strut_linear_density = self.drone.strut_beam.material_properties["mass density"] * self.drone.strut_beam.cross_section_properties["area"]
        hub_linear_density = self.drone.hub_beam.material_properties["mass density"] * self.drone.hub_beam.cross_section_properties["area"]

        self.hub_mass = hub_linear_density * (self.drone.hub_radius)
        self.arm_mass = arm_linear_density * (self.drone.nominal_rad - self.drone.hub_radius)
        self.strut_mass = strut_linear_density * self.drone.strut_length
        self.lumped_mass_frame = self.drone.blade_num * (self.arm_mass + self.strut_mass + self.hub_mass)

        self.total_mass = self.lumped_mass_frame
        self.total_mass_calculation()  # add battery masses

        # https://en.wikipedia.org/wiki/List_of_moments_of_inertia
        arm_r1 = self.drone.arm_beam.cross_section_parameters[1]/2
        arm_r2 = self.drone.arm_beam.cross_section_parameters[0]/2
        strut_r1 = self.drone.strut_beam.cross_section_parameters[1]/2
        strut_r2 = self.drone.strut_beam.cross_section_parameters[0]/2

        # moment of inertia of annuli
        self.arm_local_moment_of_inertia = np.array([[(1/2) * self.arm_mass * (arm_r1 ** 2 + arm_r2 ** 2), 0, 0],
                                                     [0, (1/12) * self.arm_mass * (3 * (arm_r1 ** 2 + arm_r2 ** 2) + self.drone.nominal_rad ** 2), 0],
                                                     [0, 0, (1/12) * self.arm_mass * (3 * (arm_r1 ** 2 + arm_r2 ** 2) + self.drone.nominal_rad ** 2)]])

        self.strut_local_moment_of_inertia = np.array([[(1/2) * self.strut_mass * (strut_r1 ** 2 + strut_r2 ** 2), 0, 0],
                                                       [0, (1/12) * self.strut_mass * (3 * (strut_r1 ** 2 + strut_r2 ** 2) + self.drone.strut_length ** 2), 0],
                                                       [0, 0, (1/12) * self.strut_mass * (3 * (strut_r1 ** 2 + strut_r2 ** 2) + self.drone.strut_length ** 2)]])

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
        self.total_moment_of_inertia = np.zeros([3, 3])
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
                      [l3, m3, n3],
                      [l2, m2, n2]])

        return t

    def calculate_transformation_matrices(self):
        for i, battery in enumerate(self.batteries):
            self.battery_transformation[i] = battery.transformation

        for i in range(self.drone.blade_num):
            theta = i * 2 * math.pi / self.drone.blade_num
            self.arm_transformation[i] = self.transformation_matrix(
                np.array([math.cos(theta), math.sin(theta), 0]),
                np.array([0, 0, 1]))

            theta_plus_1 = (i+1) * 2 * math.pi / self.drone.blade_num
            self.strut_transformation[i] = self.transformation_matrix(
                np.array([math.cos(theta_plus_1) - math.cos(theta), math.sin(theta_plus_1) - math.sin(theta), 0]),
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
            rotated_moment = self.batteries[i].local_moment_of_inertia
            displacement = self.batteries[i].position - self.center_of_mass
            translation = self.batteries[i].mass * (np.dot(displacement, displacement) * np.identity(3) - np.outer(displacement, displacement))
            # https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
            self.battery_global_moments_of_inertia[i] = rotated_moment + translation
            self.total_moment_of_inertia += self.battery_global_moments_of_inertia[i]

        for i in range(self.drone.blade_num):
            rotated_moment = self.arm_transformation[i] @ self.arm_local_moment_of_inertia @ np.transpose(self.arm_transformation[i])
            displacement = self.arm_beam_centroids[i] - self.center_of_mass
            translation = self.arm_mass * ((np.dot(displacement, displacement) * np.identity(3)) - np.outer(displacement, displacement))
            self.arm_global_moments_of_inertia[i] = rotated_moment + translation
            self.total_moment_of_inertia += self.arm_global_moments_of_inertia[i]

        for i in range(self.drone.blade_num):
            rotated_moment = self.strut_transformation[i] @ self.strut_local_moment_of_inertia @ np.transpose(self.strut_transformation[i])
            displacement = self.strut_beam_centroids[i] - self.center_of_mass
            translation = self.strut_mass * ((np.dot(displacement, displacement) * np.identity(3)) - np.outer(displacement, displacement))
            self.strut_global_moments_of_inertia[i] = rotated_moment + translation
            self.total_moment_of_inertia += self.strut_global_moments_of_inertia[i]


if __name__ == '__main__':
    from Scripts.Mass_Analysis.drone_battery import Battery
    from Scripts.drone_geometry import DroneGeometry
    from Scripts.Finite_Elements.Beam_FEA_System.beam_type import BeamType

    # 40.05 N weight
    # Mass in N-s^2/mm (Mg)
    battery_1 = Battery(4.292e-3, np.array([0, 86.75, -50]), np.array([260.5, 123.5, 63.5]))
    battery_2 = Battery(4.292e-3, np.array([0, -86.75, -50]), np.array([260.5, 123.5, 63.5]))


    # al: 0.098 lb/in^3
    # 26.6 kN/m^3 =* 10^-9 *
    arm_beam = BeamType("annulus", [29, 25], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [14, 12], "Aluminum7075-T6", "strut_beam")

    d1 = DroneGeometry("drone_test", 1000, 660, 8, arm_beam, strut_beam, [battery_1, battery_2])
    d1_mass = DroneMass(d1)

    print("\nParameters")
    print("- battery density:", battery_1.mass/(battery_1.dimensions[0] * battery_1.dimensions[1] * battery_1.dimensions[2]) * 1e12, "kg/m^3")
    print("- strut length: ", d1.strut_length, "mm")
    print("\nResults")
    print("- total mass: ", d1_mass.total_mass * 1e6, "Mg")
    print("- moment of inertia tensor: ", d1_mass.total_moment_of_inertia)
    print("- center of mass: ", d1_mass.center_of_mass, "mm")
    print("- battery ixx: ", battery_1.local_moment_of_inertia[0, 0] * 1e6, "g*mm^2")
    print("- battery iyy: ", battery_1.local_moment_of_inertia[1, 1] * 1e6, "g*mm^2")
    print("- battery izz: ", battery_1.local_moment_of_inertia[2, 2] * 1e6, "g*mm^2")

    print("- arm mass: ", d1_mass.arm_mass * 1e6, "g")
    print("- arm ixx: ", d1_mass.arm_local_moment_of_inertia[0, 0] * 1e6, "g*mm^2")
    print("- arm iyy: ", d1_mass.arm_local_moment_of_inertia[1, 1] * 1e6, "g*mm^2")
    print("- arm izz: ", d1_mass.arm_local_moment_of_inertia[2, 2] * 1e6, "g*mm^2")

    print("- strut mass: ", d1_mass.strut_mass * 1e6, "g")
    print("- strut ixx: ", d1_mass.strut_local_moment_of_inertia[0, 0] * 1e6, "g*mm^2")
    print("- strut iyy: ", d1_mass.strut_local_moment_of_inertia[1, 1] * 1e6, "g*mm^2")
    print("- strut izz: ", d1_mass.strut_local_moment_of_inertia[2, 2] * 1e6, "g*mm^2")

    print("- total mass: ", d1_mass.total_mass * 1e6, "g")
    print("- total cog z: ", d1_mass.center_of_mass, "mm")
    print("- total Ixx: ", d1_mass.total_moment_of_inertia[0, 0] * 1e6, "g*mm^2")
    print("- total Iyy: ", d1_mass.total_moment_of_inertia[1, 1] * 1e6, "g*mm^2")
    print("- total Izz: ", np.round(d1_mass.total_moment_of_inertia[2, 2] * 1e6), "g*mm^2")
    print("- total I: ",np.round(d1_mass.total_moment_of_inertia))

    print("- single beam I (on axis):\n", np.round(d1_mass.arm_global_moments_of_inertia[0] * 1e6, 3), "g*mm^2")
    print("- single beam I (on diagonal):\n", np.round(d1_mass.arm_global_moments_of_inertia[1] * 1e6, 3), "g*mm^2")
    print("- Battery I:\n", np.round((d1_mass.battery_global_moments_of_inertia[0] + d1_mass.battery_global_moments_of_inertia[1])* 1e6, 3), "g*mm^2")

    print("\nSolidWorks Results")
    print("- battery Mass:", 4292, "g")
    print("- battery ixx:", 6897422.83, "g*mm^2")
    print("- battery iyy:", 25713550.83, "g*mm^2")
    print("- battery izz:", 29726570.83, "g*mm^2")

    print("- arm beam mass:", 238.35, "g")
    print("- arm beam ixx:", 4987518.95, "g*mm^2")
    print("- arm beam iyy:", 4987518.95, "g*mm^2")
    print("- arm beam izz:", 43678.12, "g*mm^2")

    print("- strut beam mass:", 17.57, "g")
    print("- strut beam ixx:", 34675.13, "g*mm^2")
    print("- strut beam iyy:", 34675.13, "g*mm^2")
    print("- strut beam izz:", 746.60, "g*mm^2")

    print("- total mass:", 10631.36, "g")
    print("- total cog:", [0, 0, -40.37], "mm")
    print("- total I:", [164780724.15, 137813551.65, 288204656.86], "g*mm^2")

    print("- single beam I (on axis):\n", [
        [10352592.34267, 0, 2405641.75083],
        [0,  20273032.58876, 0],
        [2405641.75083, 0, 19884558.61247]])

    print("- single beam I (on diagonal):\n", [
        [10352592.34267, -9920440.24609,  1701045.59512],
        [-9920440.24609, 10352592.34267, -1701045.59512],
        [1701045.59512,  -1701045.59512, 19884558.61247]])

    print("- Battery I :\n", [
        [79190139.63839, 0, 0],
        [0,  52222967.14404, 0],
        [0, 0,  124052570.14069]])

    print("- Post 10, 20, 30 displacement on battery 1: ", [[180521391,-4235164,-354510],[-4235164,137893198,-1187895],[-354510,-11878951,304377532]])


#%%
