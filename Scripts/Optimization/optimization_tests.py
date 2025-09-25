import unittest
import numpy as np
import scipy
from Scripts.Optimization.optimization import Optimizer

from Scripts.Finite_Elements.Beam_FEA_System.beam_type import BeamType
from Scripts.Finite_Elements.Beam_FEA_System.beam_system import BeamSystem
from Scripts.Finite_Elements.material_properties import Materials

class MyTestCase(unittest.TestCase):
    def test_constraints(self):
        opt = Optimizer()
        opt.independent_properties_calc()

        opt.active_constraints = {
            "deflection": True,
            "stress": False,
            "frequency": True,
            "natural frequency": False,
            "energy": False,
            "damage": False,
        }
        '''
        
        opt.design_variables_initial_guess = {
            "arm_diameter": 50,
            "arm_thickness": 25,
            "strut_diameter": 50,
            "strut_thickness": 25,
            "strut_distance": 600,
            "hub_radius": 300,
            "hub_flange_thickness": 10,
            "hub_web_thickness": 10
        }
        '''
        res = opt.run_opt()

        x = np.array(res.x) * np.array(list(opt.design_variables_multipliers.values()))

        cons_res_opt, cons_labels_opt = opt.test_results(x)

        cons_res_int, cons_labels_int = opt.test_results(list(opt.design_variables_initial_guess.values()))

        dv_keys = list(opt.design_variables_initial_guess.keys())

        print("===== INITIAL GUESS =====")

        print("design variables")
        for i in range(len(dv_keys)):
            print(f"   {dv_keys[i]}: {list(opt.design_variables_initial_guess.values())[i]}")

        print("constraints")
        for i in range(len(cons_res_int)):
            print(f"   {cons_labels_int[i]}: {cons_res_int[i]}")

        print(f"mass w/ payload: {opt.mass(opt.design_variables_initial_guess) * 1000} kg")

        print("===== Optimized =====")
        print("design variables")
        for i in range(len(dv_keys)):
            print(f"   {dv_keys[i]}: {x[i]}")

        print("constraints")
        for i in range(len(cons_res_opt)):
            print(f"   {cons_labels_opt[i]}: {cons_res_opt[i]}")

        print(f"mass w/ payload: {res.fun * 1000} kg")

        print("OPT_OUTPUT")
        opt.test_results(x, detailed_output=True)

        print("Initial guess opt actual:")
        print(res.x)


    def test_bounds(self):
        opt = Optimizer()

        opt.independent_properties_calc()

        opt.active_constraints = {
            "deflection": False,
            "stress": False,
            "frequency": False,
            "natural frequency": False,
            "energy": False,
            "damage": False,
        }
        def bound_tst(bt):
            cons_res_int, cons_labels_int = opt.test_results(list(bt.values()))

            dv_keys = list(opt.design_variables_initial_guess.keys())


            print("design variables")
            for i in range(len(dv_keys)):
                print(f"   {dv_keys[i]}: {list(bt.values())[i]}")

            print("constraints")
            for i in range(len(cons_res_int)):
                print(f"   {cons_labels_int[i]}: {cons_res_int[i]}")

            print(f"mass: {opt.mass(bt) * 1000} kg")



        '''
        
        opt.design_variables_initial_guess = {
            "arm_diameter": 50,
            "arm_thickness": 25,
            "strut_diameter": 50,
            "strut_thickness": 25,
            "strut_distance": 600,
            "hub_radius": 300,
            "hub_flange_thickness": 10,
            "hub_web_thickness": 10
        }
        '''

        bt1 = {
            "arm_diameter": opt.boundaries["arm_diameter"][0],
            "arm_thickness": opt.boundaries["arm_thickness"][0],
            "strut_diameter": opt.boundaries["strut_diameter"][0],
            "strut_thickness": opt.boundaries["strut_thickness"][0],
            "strut_distance": opt.boundaries["strut_distance"][0],
            "hub_radius": opt.boundaries["hub_radius"][0],
            "hub_flange_thickness": opt.boundaries["hub_flange_thickness"][0],
            "hub_web_thickness": opt.boundaries["hub_web_thickness"][0]
        }
        bt2 = {
            "arm_diameter": opt.boundaries["arm_diameter"][1],
            "arm_thickness": opt.boundaries["arm_thickness"][1],
            "strut_diameter": opt.boundaries["strut_diameter"][1],
            "strut_thickness": opt.boundaries["strut_thickness"][1],
            "strut_distance": opt.boundaries["strut_distance"][1],
            "hub_radius": opt.boundaries["hub_radius"][1],
            "hub_flange_thickness": opt.boundaries["hub_flange_thickness"][1],
            "hub_web_thickness": opt.boundaries["hub_web_thickness"][1]
        }

        print("===== BOUND TEST 1 =====")
        bound_tst(bt1)

        print("===== BOUND TEST 2 =====")
        bound_tst(bt2)
    def test_optimizer(self):
        def Mass(x):
            return x[0] * x[1] * 500 * Materials["Carbon Fiber"]["mass density"] * (10 ** 6)
        def Constraints(x):
            beam = BeamType("rectangle",
                            x,
                            "Carbon Fiber",
                            "arm_beam")

            beam_system = BeamSystem("Reverse Cantilever-end load")
            beam_system.add_node(np.array([0, 0, 0]), "origin")
            beam_system.add_node(np.array([500, 0, 0]), "point_1")

            beam_system.add_beam(beam, "origin", "point_1", np.array([0, 0, 1]), "beam")

            beam_system.add_boundary_condition(0, "x", "origin")
            beam_system.add_boundary_condition(0, "y", "origin")
            beam_system.add_boundary_condition(0, "z", "origin")
            beam_system.add_boundary_condition(0, "theta x", "origin")
            beam_system.add_boundary_condition(0, "theta y", "origin")
            beam_system.add_boundary_condition(0, "theta z", "origin")

            beam_system.add_force(np.array([0,0,-40]), "point_1")

            beam_system.solve_FEA()
            beam_system.solve_static()

            disp = max(beam_system.mag_displacements)

            return disp

        consts = (
            scipy.optimize.NonlinearConstraint(Constraints, 0, 50)
        )

        results = scipy.optimize.minimize(Mass, [10,10], constraints=consts, bounds=[(5,200),(5,200)], method='SLSQP')
        print("res:", results.x)
        print("mass:", Mass(results.x))
        print("cons:", Constraints(results.x))

        print("Nathaniel Guessing")
        x = [20,5]
        print("res:", x)
        print("mass:", Mass(x))
        print("cons:", Constraints(x))

    def test_abaqus(self):
        opt = Optimizer()
        opt.simulation_settings["hub_sections"] = 1 # For simplification
        opt.independent_properties_calc()

        dv_test = opt.design_variables_initial_guess
        dv_keys = list(opt.design_variables_initial_guess.keys())

        print("===== Test Variables =====")

        print("Design Variables")
        for i in range(len(dv_keys)):
            print(f"   {dv_keys[i]}: {list(dv_test.values())[i]}")

        cons, cons_labels = opt.constraint_calculations(dv_test, opt.active_constraints_defaults, True)

        print("Abaqus Results")
        print()





if __name__ == '__main__':
    unittest.main()
