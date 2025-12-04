import unittest
import numpy as np
from Scripts.Optimization.optimization import Optimizer
from Scripts.Finite_Elements.beam_mass import BeamMass
from Scripts.Finite_Elements.drone_fea import DroneFEA
from Scripts.drone_geometry import DroneGeometry

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
        x_expanded = list(opt.reform_dvs(x).values())

        cons_res_opt, cons_labels_opt = opt.test_results(x)

        dv_keys = list(opt.design_variables_initial_guess.keys())

        print("===== Optimized =====")
        print("design variables")
        for i in range(len(dv_keys)):
            print(f"   {dv_keys[i]}: {x_expanded[i]}")

        print("constraints")
        for i in range(len(cons_res_opt)):
            print(f"   {cons_labels_opt[i]}: {cons_res_opt[i]}")

        print(f"mass w/ payload: {res.fun * 1000} kg")

        print("OPT_OUTPUT")
        opt.test_results(x, detailed_output=True)

        print("Initial guess opt actual:")
        print(res.x)

        geo = opt.drone_geometry
        sec = 50
        geo_update = DroneGeometry( geo.name + "_update", geo.nominal_rad, geo.hub_radius, geo.strut_pos, geo.blade_num, geo.arm_beam, geo.strut_beam, geo.hub_beam, sec, geo.batteries)
        geo_update.calc_node_coords()

        d1_fea = DroneFEA(geo_update)
        d1_fea.create_drone_nodes()
        d1_fea.create_drone_beams()
        d1_fea.boundary_conditions()
        print("beam witdhs (mm): ", geo_update.beam_widths)


        mass_analysis = BeamMass(d1_fea.beam_system)
        print("Beam Count: ", len(d1_fea.beam_system.beams))
        print("Mass (kg): ", mass_analysis.total_mass * 1000)
        print("Moment of Inertia (kg/mm^2):", mass_analysis.total_moment_of_inertia * 1000)
        print("Center of Mass (mm): ", mass_analysis.center_of_mass)
        # print(d1_fea.beam_system.beam_names)
        arm_mass =  mass_analysis.beam_mass[d1_fea.beam_system.select_element('outer_arm_beam_0')] + mass_analysis.beam_mass[d1_fea.beam_system.select_element('inner_arm_beam_0')]
        print("Arm Mass (kg): ", arm_mass  * 1000)
        print("Strut Mass (kg): ", mass_analysis.beam_mass[d1_fea.beam_system.select_element('strut_beam_0')] * 1000)
        hub_mass = 0
        for i in range(8):
            for j in range(sec):
                hub_mass += mass_analysis.beam_mass[d1_fea.beam_system.select_element(f'hub_beam_{j}_{i}')]

        print("Hub Mass (kg): ", hub_mass * 1000)
        print("hub segment length (mm):",  d1_fea.beam_system.beams[d1_fea.beam_system.select_element(f'hub_beam_0_0')].length)
        print("hub segment length * seg (mm):",  d1_fea.beam_system.beams[d1_fea.beam_system.select_element(f'hub_beam_0_0')].length * sec)
        print("hub segment 0 area (mm^2):",  d1_fea.beam_system.beams[d1_fea.beam_system.select_element(f'hub_beam_0_0')].beam_type.cross_section_properties["area"])

if __name__ == '__main__':
    unittest.main()
