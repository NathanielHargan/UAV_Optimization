from Scripts.Finite_Elements.cross_section_properties import cross_section_circle, cross_section_annulus, \
    cross_section_rectangle, cross_section_hexagon, cross_section_i_beam
from Scripts.Finite_Elements.material_properties import Materials
from Scripts.Mass_Analysis import mass_properties as MassProperties


class BeamType:
    def __init__(self, cross_section, cross_section_parameters, material, name=None):
        self.cross_section = cross_section
        self.cross_section_parameters = cross_section_parameters
        self.material_properties = Materials[material]

        # Name is optional so
        if name is None:
            self.name = "Beam_Type_" + cross_section + "_" + material
        else:
            self.name = name

        self.cross_section_properties = self.cross_section_properties_init()


        # for shape functions
        eik_z = (self.material_properties["elastic modulus"] *
                 self.cross_section_properties["second moment of area z"] *
                 self.cross_section_properties["transverse shear deflection constant z"])

        eik_y = (self.material_properties["elastic modulus"] *
                 self.cross_section_properties["second moment of area y"] *
                 self.cross_section_properties["transverse shear deflection constant y"])

        ga = self.material_properties["shear modulus"] * self.cross_section_properties["area"]

        self.g_y = eik_y/ga
        self.g_z = eik_z/ga



    def mass_moment(self, L):
        if self.cross_section.lower() == "circle":
            return MassProperties.mass_moment_circle(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "annulus":
            return MassProperties.mass_moment_annulus(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "rectangle":
            return MassProperties.mass_moment_rectangle(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "hexagon":
            return MassProperties.mass_moment_hexagon(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        elif self.cross_section.lower() == "i_beam":
            return MassProperties.mass_moment_i_beam(
                self.cross_section_parameters,
                L,
                self.material_properties["mass density"])
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'


    def cross_section_properties_init(self):
        if self.cross_section.lower() == "circle":
            return cross_section_circle(self.cross_section_parameters)
        elif self.cross_section.lower() == "annulus":
            return cross_section_annulus(self.cross_section_parameters)
        elif self.cross_section.lower() == "rectangle":
            return cross_section_rectangle(self.cross_section_parameters)
        elif self.cross_section.lower() == "hexagon":
            return cross_section_hexagon(self.cross_section_parameters)
        elif self.cross_section.lower() == "i_beam":
            return cross_section_i_beam(self.cross_section_parameters)
        else:
            print('incorrect cross section in beam ' + self.name)
            return 'error'
