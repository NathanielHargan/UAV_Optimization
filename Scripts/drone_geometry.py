import numpy as np
import math


class DroneGeometry:
    def __init__(self, name, nominal_rad, hub_radius, strut_pos, blade_num, arm_beam, strut_beam, hub_beam, hub_sections, batteries=[]):
        self.name = name  # Name of the drone

        self.nominal_rad = nominal_rad
        self.hub_radius = hub_radius
        self.strut_pos = strut_pos
        self.blade_num = blade_num

        self.hub_beam = hub_beam
        self.arm_beam = arm_beam
        self.strut_beam = strut_beam

        self.hub_sections = hub_sections

        self.hub_edge_length = 2 * hub_radius * math.sin(math.pi/blade_num)

        theta = 2 * math.pi / blade_num
        x_next = math.cos(theta) * strut_pos
        y_next = math.sin(theta) * strut_pos
        x = (x_next + strut_pos)/2
        y = y_next/2

        self.strut_end_x_pos = x
        self.strut_end_y_pos = y
        self.hub_nodes_coords = np.zeros([blade_num, hub_sections+1, 3])

        self.beam_widths = np.zeros([hub_sections])
        h = math.cos(math.pi/blade_num) * hub_radius
        bh_ratio = self.hub_edge_length / h
        h_seg = h/hub_sections
        for i in range(hub_sections):
            h_i = h * (i/hub_sections)
            h_ip1 = h * ((i+1)/hub_sections)
            w_i = h_i * bh_ratio
            w_ip1 = h_ip1 * bh_ratio
            # dist from h_i+1
            centroid_dist = (h_seg /3) * (w_ip1 + 2 * w_i) / (w_ip1 + w_i)
            width_at_centroid = bh_ratio * (h_ip1 - centroid_dist)
            self.beam_widths[i] = width_at_centroid



        self.outer_nodes_coords = np.zeros([blade_num,3])
        self.strut_nodes_coords = np.zeros([blade_num,3])
        self.strut_centroid_coords = np.zeros([blade_num,3])
        self.outer_centroid_coords = np.zeros([blade_num,3])
        self.calc_node_coords()

        self.strut_length = np.linalg.norm(self.strut_nodes_coords[0] - self.strut_nodes_coords[1])

        self.arm_length = nominal_rad - hub_radius

        self.batteries = batteries

        self.projected_surface_area = 0
        self.calc_projected_surface_area()


    def calc_node_coords(self):
        for i in range(self.blade_num):
            theta = i * (2 * math.pi / self.blade_num)
            theta_next = (i+1) * (2 * math.pi / self.blade_num)
            self.outer_nodes_coords[i, 0] = math.cos(theta) * self.nominal_rad
            self.outer_nodes_coords[i, 1] = math.sin(theta) * self.nominal_rad
            self.strut_nodes_coords[i, 0] = math.cos(theta) * self.strut_pos
            self.strut_nodes_coords[i, 1] = math.sin(theta) * self.strut_pos

            r_list = np.linspace(0,  self.hub_radius, num=self.hub_sections+1, endpoint=True)

            for j, r in enumerate(r_list):
                self.hub_nodes_coords[i, j, 0] = math.cos(theta) * r
                self.hub_nodes_coords[i, j, 1] = math.sin(theta) * r

            self.outer_centroid_coords[i, 0] = math.cos(theta) * self.nominal_rad/2
            self.outer_centroid_coords[i, 1] = math.sin(theta) * self.nominal_rad/2

            self.strut_centroid_coords[i, 0] = (math.cos(theta) + math.cos(theta_next)) * self.strut_pos / 2
            self.strut_centroid_coords[i, 1] = (math.sin(theta) + math.sin(theta_next)) * self.strut_pos / 2

    def calc_projected_surface_area(self):
        battery_sa = 0
        for battery in self.batteries:
            battery_sa += battery.dimensions[0] * battery.dimensions[1]

        strut_sa = self.blade_num * (self.strut_beam.cross_section_parameters[0]) * self.strut_length
        arm_sa = self.blade_num * (self.arm_beam.cross_section_parameters[0]) * (self.nominal_rad - self.hub_radius)
        hub_sa = self.blade_num * self.hub_edge_length * self.hub_radius * math.cos(math.pi/self.blade_num)

        total_sa = arm_sa + strut_sa + battery_sa + hub_sa
        self.projected_surface_area = total_sa













    #%%

#%%
