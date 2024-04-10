import numpy as np
from Scripts.cross_section_properties import cross_section_circle
from Scripts.cross_section_properties import cross_section_annulus
from Scripts.material_properties import Materials

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
            name = "Node_" + str(len(self.nodes))

        self.node_names.append(name)

    def add_beam(self, beam_type, start_node, end_node, name=None):
        if type(start_node) is int:
            start_node_coords = self.nodes[start_node]
            end_node_coords = self.nodes[end_node]
        elif type(start_node) is str:
            # Finds node with the same name as the inputted string
            start_node_coords = self.nodes[self.node_names.index(start_node)]
            end_node_coords = self.nodes[self.node_names.index(end_node)]
        else:
            start_node_coords = [0, 0, 0]
            end_node_coords = [0, 0, 0]
            print("start or end nodes must be string or integer")

        if name is None:
            name = "Beam_" + str(len(self.beams))

        new_beam = Beam(beam_type, start_node_coords, end_node_coords, name)
        self.beams.append(new_beam)



class Beam:
    def __init__(self, beam_type, start_node_coords, end_node_coords, name=None):
        self.name = name  # String with the beam name
        self.beam_type = beam_type

        # Coordinates of the nodes
        self.start_node_coords = start_node_coords
        self.end_node_coords = end_node_coords

        self.length = np.linalg.norm(end_node_coords - start_node_coords)


class BeamType:
    def __init__(self, cross_section, cross_section_parameters, material, name=None):
        self.cross_section = cross_section
        self.cross_section_parameters = cross_section_parameters
        self.material = Materials[material]

        # Name is optional so
        if name is None:
            self.name = "Beam_Type_" + cross_section + "_" + material
        else:
            self.name = name

        self.cross_section_properties = self.cross_section_properties_init()

    def cross_section_properties_init(self):
        if self.cross_section.lower() == "circle":
            return cross_section_circle(self.cross_section_parameters)
        elif self.cross_section.lower() == "annulus":
            return cross_section_annulus(self.cross_section_parameters)
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'

#%%
