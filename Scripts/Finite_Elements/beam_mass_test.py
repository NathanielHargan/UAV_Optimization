import unittest
import numpy as np
from Scripts.Optimization.optimization import Optimizer
from Scripts.Finite_Elements.beam_mass import BeamMass
from Scripts.Finite_Elements.drone_fea import DroneFEA

class MyTestCase(unittest.TestCase):
    def test_mass_moment_opt(self):
        opt = Optimizer()
        opt.independent_properties_calc()
        opt.active_constraints = {
            "deflection": True,
            "stress": True,
            "frequency": True,
            "natural frequency": False,
            "energy": False,
            "damage": False,
        }

        opt.design_variables_initial_guess = {
            "arm_diameter": 50.0,
            "arm_thickness": 0.6,
            "strut_diameter": 50.0,
            "strut_thickness": 0.6,
            "strut_distance": 475,
            "hub_radius": 250,
            "hub_flange_thickness": 0.7,
            "hub_web_thickness": 0.7
        }

        res = opt.run_opt()

        x = np.array(res.x)
        print(res.x)
        opt.test_results(x, detailed_output=True)

        geo = opt.drone_geometry


        d1_fea = DroneFEA(geo)
        d1_fea.create_drone_nodes()
        d1_fea.create_drone_beams()
        d1_fea.boundary_conditions()

        mass_analysis = BeamMass(d1_fea.beam_system)
        print(mass_analysis.total_mass)
        print(mass_analysis.total_moment_of_inertia)
        print(mass_analysis.center_of_mass)

        print(mass_analysis.beam_moment_of_inertia)
        print(mass_analysis.beam_mass)
        print(mass_analysis.beam_centroids)



if __name__ == '__main__':
    unittest.main()
