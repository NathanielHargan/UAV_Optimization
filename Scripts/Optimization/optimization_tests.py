import unittest
from Scripts.Optimization.optimization import Optimizer

class MyTestCase(unittest.TestCase):
    def test_constraints(self):
        opt = Optimizer()
        opt.independent_properties_calc()

        opt.active_constraints = {
            "deflection": True,
            "stress": True,
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

        cons_res_opt, cons_labels_opt = opt.test_results(res.x)

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
            print(f"   {dv_keys[i]}: {res.x[i]}")

        print("constraints")
        for i in range(len(cons_res_opt)):
            print(f"   {cons_labels_opt[i]}: {cons_res_opt[i]}")

        print(f"mass w/ payload: {res.fun * 1000} kg")






if __name__ == '__main__':
    unittest.main()
