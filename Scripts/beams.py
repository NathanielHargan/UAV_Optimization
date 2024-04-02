import numpy as np
from Scripts import cross_section_properties as csp


class BeamSystem:
    def __init__(self, name):
        self.name = name  # A string with the name of the system
        self.nodes = np.empty((0, 3))  # List of nodes objects
        self.node_names = []  # List of nodes names
        self.beams = []  # List of beam objects

    def add_node(self, x, y, z, node_name = None):
        coords = np.array([[x, y, z]])
        self.nodes = np.concatenate((self.nodes, coords), axis=0)
        if node_name is None:
            node_name = "Node_" + str(len(self.nodes))

        self.node_names.append(name)

    def add_beam(self, start_node, end_node, cross_section, cross_section_parameters, material, beam_name=None):
        start_node_coords = self.nodes[start_node]
        end_node_coords = self.nodes[end_node]

        if beam_name is None:
            beam_name = "Beam_" + str(len(self.beams))

        self.beams.append(
            Beam(beam_name,
                 start_node_coords,
                 end_node_coords,
                 cross_section,
                 cross_section_parameters,
                 material))


class Beam:
    def __init__(self, name, start_node_coords, end_node_coords, cross_section, cross_section_parameters, material):
        self.name = name  # String with the beam name
        self.cross_section = cross_section
        self.cross_section_parameters = cross_section_parameters
        self.material = material
        self.start_node_coords = start_node_coords
        self.end_node_coords = end_node_coords
        self.cross_section_properties = {}
        self.cross_section_properties()

    def cross_section_properties(self):
        if self.cross_section.lower() == "circle":
            self.cross_section_properties = csp.cross_section_circle(self.cross_section_parameters)
        if self.cross_section.lower() == "annulus":
            self.cross_section_properties = csp.cross_section_annulus(self.cross_section_parameters)

