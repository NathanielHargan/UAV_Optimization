import numpy as np
import math


class Drone:
    def __init__(self, name, mass_drone_body):
        self.name = name # Name of the drone
        self.mass_drone_body = mass_drone_body # Mass of drone without payload

        self.nodes = np.empty((0, 3))
        self.beams = np.empty((0, 2), int)

    def add_node(self, x, y, z):
        coords = np.array([[x, y, z]])
        self.nodes = np.concatenate((self.nodes, coords), axis=0)

    def add_beam(self, a, b):
        v = np.array([[a, b]])
        self.beams = np.concatenate((self.beams, v), axis=0)

    def create_drone(self, hub_rad_out, nominal_rad, blade_num):
        blade_1_theta = 2*math.pi/blade_num

        # Midpoint between this node and next node
        # hub_rad_out * cos(blade_1_theta) - hub_rad_out * cos(0)
        x_calc = (hub_rad_out*math.cos(blade_1_theta) + hub_rad_out)/2

        # hub_rad_out * sin(blade_1_theta) - hub_rad_out * sin(0)
        y_calc = (hub_rad_out*math.sin(blade_1_theta))/2
        inner_hub_radius = np.linalg.norm([x_calc,y_calc])

        for theta in np.linspace(0, 2*math.pi, blade_num+1):
            x_out = math.cos(theta) * hub_rad_out
            y_out = math.sin(theta) * hub_rad_out
            z_out = 0
            self.add_node(x_out, y_out, z_out)  # outer hub node

            # Midpoint between this node and next node
            x_in = math.cos(theta + math.pi/blade_num) * inner_hub_radius
            y_in = math.sin(theta + math.pi/blade_num) * inner_hub_radius
            z_in = 0

            self.add_node(x_in, y_in, z_in)  # inner hub node

            x_ex = math.cos(theta + math.pi/blade_num) * nominal_rad
            y_ex = math.sin(theta + math.pi/blade_num) * nominal_rad
            z_ex = 0

            self.add_node(x_ex, y_ex, z_ex)  # external node

        for i in range(0, len(self.nodes)-3, 3):
            outer_hub_node_index = i
            inner_hub_node_index = i+1
            external_hub_node_index = i+2

            self.add_beam(outer_hub_node_index,inner_hub_node_index)
            self.add_beam(inner_hub_node_index,external_hub_node_index)
            self.add_beam(inner_hub_node_index,outer_hub_node_index+3)

        self.add_beam(len(self.nodes)-2,0)








    

        

#%%
