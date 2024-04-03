import numpy as np
from Scripts.cross_section_properties import cross_section_circle
from Scripts.cross_section_properties import cross_section_annulus


class BeamSystem:
    def __init__(self, name):
        self.system_name = name  # A string with the name of the system
        self.nodes = np.empty((0, 3))  # List of nodes objects
        self.node_names = []  # List of nodes names
        self.beams = []  # List of beam objects

    def add_node(self, x, y, z, name=None):
        coords = np.array([[x, y, z]])
        self.nodes = np.concatenate((self.nodes, coords), axis=0)

        if name is None:
            node_name = "Node_" + str(len(self.nodes))

        self.node_names.append(name)

    def add_beam(self, start_node, end_node, cross_section, cross_section_parameters, material, name=None):
        start_node_coords = [0, 0, 0]
        end_node_coords = [0, 0, 0]

        if type(start_node) is int:
            start_node_coords = self.nodes[start_node]
            end_node_coords = self.nodes[end_node]
        elif type(start_node) is str:
            # Finds node with the same name as the inputted string
            start_node_coords = self.nodes[self.node_names.index(start_node)]
            end_node_coords = self.nodes[self.node_names.index(end_node)]
        else:
            print("start or end nodes must be string or integer")

        if name is None:
            name = "Beam_" + str(len(self.beams))

        self.beams.append(
            Beam(name,
                 start_node_coords,
                 end_node_coords,
                 cross_section,
                 cross_section_parameters,
                 material))


class Beam:
    def __init__(self, name, start_node_coords, end_node_coords, cross_section, cross_section_parameters, material):
        self.name = name  # String with the beam name
        self.cross_section = cross_section  # Type of cross-section "Circle" or "Annulus"

        # [r1] for circle [r1,r2] for annulus
        self.cross_section_parameters = cross_section_parameters

        # Will add this soon. I do not know how I will structure this data. I will figure it out.
        self.material = material

        # Coordinates of the nodes
        self.start_node_coords = start_node_coords
        self.end_node_coords = end_node_coords

        # Empty dictionary which is filled with the
        self.cross_section_properties = {}
        self.cross_section_properties = self.cross_section_properties_init()

    def cross_section_properties_init(self):
        if self.cross_section.lower() == "circle":
            return cross_section_circle(self.cross_section_parameters)
        elif self.cross_section.lower() == "annulus":
            return cross_section_annulus(self.cross_section_parameters)
        else:
            print('incorrect cross section in beam ' + self.name)