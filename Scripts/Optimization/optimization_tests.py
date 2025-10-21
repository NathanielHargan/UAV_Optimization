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
            "stress": True,
            "frequency": True,
            "natural frequency": False,
            "energy": False,
            "damage": False,
        }
        opt.design_variables_initial_guess["strut_distance"] = 575
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

        x = np.array(res.x)

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
        '''
        NO MULT
        design variables
           arm_diameter: 49.99999999999997
           arm_thickness: 0.3211725866649314
           strut_diameter: 50.0
           strut_thickness: 0.3591731107507465
           strut_distance: 495.21354488405416
           hub_radius: 100.0
           hub_flange_thickness: 0.7000000000000001
           hub_web_thickness: 0.7000000000000022
        constraints
           deflection: 0.9999999999985973
           stress: 0.9568074766549263
           frequency: 0.27331598356573267
           energy: 0.3260688812472388
           damage: 0.0
        '''
        # mass: 0.8275629847103303

        '''
        Starting 600
        ===== Optimized =====
            design variables
               arm_diameter: 50.0
               arm_thickness: 0.29458154397114766
               strut_diameter: 49.99999999999961
               strut_thickness: 0.3462775923814778
               strut_distance: 598.6565815937422
               hub_radius: 104.39581041184776
               hub_flange_thickness: 0.7
               hub_web_thickness: 0.7
            constraints
               deflection: 1.000000497957219
               stress: 1.0159101537875042
               frequency: 1.000193766034515
               energy: 0.3260837080288574
               damage: 0.0
            mass w/o payload: 0.8272803134225961 kg
        '''

        '''
            RAISED ITER LIMIT
            ===== Optimized =====
            design variables
               arm_diameter: 50.0
               arm_thickness: 0.2963381407683932
               strut_diameter: 49.999999999999865
               strut_thickness: 0.34910466145795416
               strut_distance: 601.8658063700772
               hub_radius: 100.0
               hub_flange_thickness: 0.7000000000000013
               hub_web_thickness: 0.700000000000012
            constraints
               deflection: 1.0000000000000007
               stress: 1.0131786520729813
               frequency: 0.9999999999995449
               energy: 0.32606836206828305
               damage: 0.0
            mass w/ payload: 0.8265502599010393 kg
            OPT_OUTPUT
            === Optimizer Results ===
            Wingspan: 1127.3025260953375 mm
            Force: -55.162406246625004 N
            Hub Beam Width: 55.22847498307934 mm
            Tip Displacement: 50.000000000000036 mm
            Strut Center Displacement: 15.242743721339828 mm
            Stress Strut: 37.385234634421295 MPa
            Stress Hub: 62.393099678285 MPa
            Stress Outer Beam: 50.7078967987886 MPa
            Tip Displacement Amplitude: 74.99999999996592 mm
            Propeller Frequency 355.51173598844576 rad/s
            Initial guess opt actual:
            [5.00000000e+01 2.96338141e-01 5.00000000e+01 3.49104661e-01
             6.01865806e+02 1.00000000e+02 7.00000000e-01 7.00000000e-01]
        '''

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






if __name__ == '__main__':
    unittest.main()
