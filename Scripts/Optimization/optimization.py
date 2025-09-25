
import pandas as pd
import pylife.stress.rainflow as RF
from pylife.strength import fatigue
import math
import numpy as np
import scipy
import Scripts.unit_conversions as uc
from Scripts.Finite_Elements.material_properties import Materials
from Scripts.Finite_Elements.Beam_FEA_System.beam_type import BeamType
from Scripts.Finite_Elements.material_properties import Gravity
from Scripts.Mass_Analysis.drone_battery import Battery
from Scripts.Mission_Profile.drone_power import BatterySystem
from Scripts.drone_geometry import DroneGeometry
from Scripts.Finite_Elements.drone_fea import DroneFEA
from Scripts.Mass_Analysis.drone_mass import DroneMass
from Scripts.Mission_Profile.drone_power import DronePower
from Scripts.Mission_Profile.mission_profile import MissionProfileLinearAcc

class Optimizer:
    """
        Optimizer class

        Given a set of constraints, design variables, this class will calculate the optimal design variables for the drone.

        Attributes:
        active_constraints: dict
            dictionary of active constraints
            default: {
                "deflection": True,
                "stress": True,
                "natural frequency": True,
                "energy": True,
                "damage": True,
            }

        simulation_settings: dict
            dictionary of simulation settings
            default: {
                "hub_sections":5,
                "transient_timesteps": 1000,
                "energy_timesteps": 1000
            }

        parameters: dict
            defaults = {
                "blade_num":8,
                "propeller_radius":406.4,
                "propeller_spacing":25,
                "payload":45 * 10 ** -3,
                "air_density":1.225e-12,
                "drag_coeff":0.256,
                "propeller_max_RPM": 3068,
                "battery_mass": 4.292e-3,
                "battery_dimensions": np.array([260.5, 123.5, 63.5]),
                "battery_volts": 22.8,
                "battery_mah": 40000,
                "rev_up_time": 5,
                "rev_down_time": 5,
                "gravity": Gravity,
            }


        Methods
        run_opt():
            returns the scipy optimization minimization for a given set of constraints

        test_results(results):
            given the results of run_opt(), this function will check every constraint and return their values

    """
    def __init__(self):
        self.active_constraints_defaults = {
            "deflection":True,
            "stress":True,
            "natural frequency":False,
            "energy":True,
            "damage":True,
            "frequency":True
        }

        self.simulation_settings_defaults = {
            "hub_sections":5,
            "transient_timesteps": 1000,
            "energy_timesteps": 1000,
            "mass_multiplier": 1000
        }

        self.parameters_defaults = {
            "blade_num":8,
            "propeller_radius":406.4,
            "propeller_spacing":25,
            "payload":45 * 10 ** -3,
            "air_density":1.225e-12,
            "drag_coeff":0.256,
            "propeller_max_RPM": 3068,
            "battery_mass": 4.292e-3,
            "battery_dimensions": np.array([260.5, 123.5, 63.5]),
            "battery_volts": 22.8,
            "battery_mah": 40000,
            "rev_up_time": 5,
            "rev_down_time": 5,
            "gravity": Gravity,
            "hub_height": 64, # https://genstattu.com/tattu-40000mah-6s-10c-22-8v-high-voltage-uav-lipo-battery-pack-with-as150-as150-plug/?srsltid=AfmBOoqN5rYTguPKErrR2RGJXuesTt9UqwIoJ9aieuGaF-ipxM2Kz1VW
        }

        self.constraints_constants_defaults = {
            "allowable_deflection": 50,
            "stress_FOS": 4,
            "allowable_natural_frequency_FOS":1.2,
            "allowable_amplitude_to_disp_ratio": 1.5,
            "energy_remaining":0.2,
            "mission_count_damage":10000,
        }

        self.boundaries_defaults = {
            "arm_diameter": (9, 50),
            "arm_thickness": (0.254, 3.175),
            "strut_diameter": (9, 50),
            "strut_thickness": (0.254, 3.175),
            "strut_distance": (400, 1500),
            "hub_radius": (100, 400),
            "hub_flange_thickness": (0.7, 12.7),
            "hub_web_thickness": (0.7, 3.175)
        }

        self.design_variables_multipliers = {
            "arm_diameter": 100,
            "arm_thickness": 1,
            "strut_diameter": 100,
            "strut_thickness": 1,
            "strut_distance": 1000,
            "hub_radius": 100,
            "hub_flange_thickness": 1,
            "hub_web_thickness": 1
        }

        # Normalize these
        self.design_variables_initial_guess_defaults = {
            "arm_diameter": 29,
            "arm_thickness": 2,
            "strut_diameter": 29,
            "strut_thickness": 2,
            "strut_distance": 500,
            "hub_radius": 200,
            "hub_flange_thickness": 2,
            "hub_web_thickness": 2
        }


        self.active_constraints = self.parameters_defaults.copy()
        self.parameters = self.parameters_defaults.copy()
        self.design_variables_initial_guess = self.design_variables_initial_guess_defaults.copy()
        self.constraints_constants = self.constraints_constants_defaults.copy()
        self.simulation_settings = self.simulation_settings_defaults.copy()
        self.boundaries = self.boundaries_defaults.copy()



    def independent_properties_calc(self):

        self.mission_profile = self.mission_profile_calc()

        self.total_mission_duration = self.mission_profile.duration + self.parameters["rev_up_time"] + self.parameters["rev_down_time"]

        self.battery_system = BatterySystem("Tattu 40000mAh 6S 10C 22.8V High Voltage UAV Lipo Battery Pack with AS150+AS150",
            self.parameters["battery_mah"],
            self.parameters["battery_volts"],
            2)

        self.energy_calc_timesteps = np.linspace(0, self.total_mission_duration, self.simulation_settings["energy_timesteps"])

        self.completely_reversed_stress_amplitude = 1466 * (2 * 1e6) ** (-0.143)

        self.woehler_curve_data_carbon_fiber = pd.Series({
                'SD': self.completely_reversed_stress_amplitude,
                'ND': 1e6,
                'k_1': -0.143
        }) # Data soon

        self.mass_battery_1 = Battery(self.parameters["battery_mass"],
                                      np.array([0, 86.75, -50]),
                                      self.parameters["battery_dimensions"])

        self.mass_battery_2 = Battery(self.parameters["battery_mass"],
                                      np.array([0, -86.75, -50]),
                                      self.parameters["battery_dimensions"])



    def mission_profile_calc(self):
        ''' calculates the mission profile for the drone '''
        m1 = MissionProfileLinearAcc("profile")

        m1.add_segment(15, "Takeoff")
        m1.add_segment(10, "Climb1")
        m1.add_segment(10, "Climb2")
        m1.add_segment(25, "cruise1")
        m1.add_segment(25, "cruise2")
        m1.add_segment(10, "decent1")
        m1.add_segment(10, "decent2")
        m1.add_segment(15, "land")

        # Takeoff
        m1.add_constraint(0, "x", 0)
        m1.add_constraint(0, "y", 0)
        m1.add_constraint(0, "vx", 0)
        m1.add_constraint(0, "vy", 0)
        m1.add_constraint(0, "ax", 0)
        m1.add_constraint(0, "ay", 0)

        # Cruise
        m1.add_constraint(550000, "x", 3)
        m1.add_constraint(100000, "y", 3)
        m1.add_constraint(0, "ax", 3)
        m1.add_constraint(0, "vy", 3)
        m1.add_constraint(1250000, "x", 4)
        m1.add_constraint(105000, "y", 4)
        m1.add_constraint(2000000, "x", 5)
        m1.add_constraint(0, "ax", 5)
        m1.add_constraint(110000, "y", 5)
        m1.add_constraint(0, "vy", 5)

        # Decent
        m1.add_constraint(2400000, "x", 8)
        m1.add_constraint(0, "y", 8)
        m1.add_constraint(0, "vx", 8)
        m1.add_constraint(0, "vy", 8)
        m1.add_constraint(0, "ax", 8)
        m1.add_constraint(0, "ay", 8)

        m1.solve_kinematics()

        return m1



    def force_transient(self, t, m, sa):
        '''
            forces applied to the drone to follow the mission profile.
            args:
                t: time (s)
                sa: surface area (mm^2)
                m: mass (Mg)

            returns:
                a numpy array describing force [x, y, z, mxx, myy, mzz] (N / N*mm)
        '''
        # sa = surface_area
        lift_thrust = m / self.parameters["blade_num"]
        rev_time = self.parameters["rev_up_time"]
        if t < rev_time:
            return np.array([0, 0, (t / rev_time) * lift_thrust, 0, 0, 0])
        elif t < rev_time + self.mission_profile.duration:
            # Calculate force per propeller
            p, v, a = self.mission_profile.time_solve(t - rev_time)
            return np.array([0, 0, np.sign(v[1]) * (v[1] ** 2) * sa * self.parameters["air_density"] * self.parameters["drag_coeff"] + (
                        a[1] + self.parameters["gravity"]) * m / self.parameters["blade_num"], 0, 0, 0])
        else:
            return np.array([0, 0, lift_thrust - ((t - rev_time - self.mission_profile.duration) / rev_time) * (lift_thrust), 0, 0, 0])

    def force_oscillation_transient(self, t, dpm, m, sa):
        '''
            Returns a oscillating force describing the vibration from the propeller.
            args:
                t: time (s)
                dpm: drone power module (see drone_power_module.py)

            returns:
                a scalar force (N)

        '''
        return math.sin(dpm.freq_calc(self.parameters["propeller_max_RPM"], 16, t) * t) * self.force_transient(t, m, sa)[2] * 0.01

    def total_force_y(self, t, m, sa):
        return self.force_transient(t, m, sa)[2] * 8

    def mass(self, design_variables):
        wingspan = (self.parameters["propeller_spacing"] + self.parameters["propeller_radius"]) / math.sin(
            math.pi / self.parameters["blade_num"])

        arm_inner_diameter = design_variables["arm_diameter"] - design_variables["arm_thickness"] * 2
        strut_inner_diameter = design_variables["strut_diameter"] - design_variables["strut_thickness"] * 2
        hub_side_length = 2 * design_variables["hub_radius"] * math.tan(math.pi / self.parameters["blade_num"])

        arm_beam_properties = BeamType("annulus",
                                       [design_variables["arm_diameter"], arm_inner_diameter],
                                       "Carbon Fiber",
                                       "arm_beam")

        strut_beam_properties = BeamType("annulus",
                                         [design_variables["strut_diameter"],
                                          strut_inner_diameter],
                                         "Carbon Fiber",
                                         "strut_beam")

        hub_beam_properties = BeamType("i_beam",
                                       [hub_side_length * (2 / 3),
                                        self.parameters["hub_height"] + 2 * design_variables["hub_flange_thickness"],
                                        design_variables["hub_flange_thickness"],
                                        design_variables["hub_web_thickness"]],
                                       "Carbon Fiber",
                                       "hub_beam")

        drone_geometry = DroneGeometry("drone",
                                       wingspan,
                                       design_variables["hub_radius"],
                                       design_variables["strut_distance"],
                                       self.parameters["blade_num"],
                                       arm_beam_properties,
                                       strut_beam_properties,
                                       hub_beam_properties,
                                       self.simulation_settings["hub_sections"],
                                       [self.mass_battery_1, self.mass_battery_2])

        drone_mass = DroneMass(drone_geometry).lumped_mass_frame

        return drone_mass
    def constraint_calculations(self, design_variables, active_constraints, detailed_output=False):
        constraints = []
        constraint_labels = []
        # arm length calculation
        wingspan = (self.parameters["propeller_spacing"] + self.parameters["propeller_radius"]) / math.sin(math.pi/self.parameters["blade_num"])

        arm_inner_diameter = design_variables["arm_diameter"] - design_variables["arm_thickness"] * 2
        strut_inner_diameter = design_variables["strut_diameter"] - design_variables["strut_thickness"] * 2
        hub_side_length = 2 * design_variables["hub_radius"] * math.tan(math.pi/self.parameters["blade_num"])

        arm_beam_properties = BeamType("annulus",
                                       np.array([design_variables["arm_diameter"], arm_inner_diameter]),
                                       "Carbon Fiber",
                                       "arm_beam")

        strut_beam_properties = BeamType("annulus",
                                         np.array([design_variables["strut_diameter"],
                                          strut_inner_diameter]),
                                         "Carbon Fiber",
                                         "strut_beam")

        hub_beam_properties = BeamType("i_beam",
                                       np.array([hub_side_length,
                                        self.parameters["hub_height"] + 2 * design_variables["hub_flange_thickness"],
                                        design_variables["hub_flange_thickness"],
                                        design_variables["hub_web_thickness"]]),
                                       "Carbon Fiber",
                                       "hub_beam")

        drone_geometry = DroneGeometry("drone",
                                       wingspan,
                                       design_variables["hub_radius"],
                                       design_variables["strut_distance"],
                                       self.parameters["blade_num"],
                                       arm_beam_properties,
                                       strut_beam_properties,
                                       hub_beam_properties,
                                       self.simulation_settings["hub_sections"],
                                       [self.mass_battery_1, self.mass_battery_2])

        mass_module = DroneMass(drone_geometry)
        drone_mass = mass_module.total_mass + self.parameters["payload"]

        # FEA
        d1_fea = DroneFEA(drone_geometry)
        d1_fea.create_drone_slice_nodes()
        d1_fea.create_drone_slice_beams()
        d1_fea.boundary_conditions_slice()

        static_force = -self.parameters["gravity"] * (self.parameters["payload"] / self.parameters["blade_num"])
        d1_fea.beam_system.add_force(np.array([0, 0, static_force]), "outer_node")
        d1_fea.solve_fea()
        d1_fea.solve_static()

        if active_constraints["deflection"] or active_constraints["stress"] or active_constraints["frequency"]:
            d1_fea.solve_static()

        if active_constraints["deflection"]:
            max_disp = max(d1_fea.beam_system.mag_displacements)
            FOS_disp =  self.constraints_constants["allowable_deflection"] / max_disp
            deflection_constraint = max_disp / self.constraints_constants["allowable_deflection"]
            constraints.append(deflection_constraint)
            constraint_labels.append("deflection")

        if active_constraints["stress"]:
            d1_fea.solve_failure()
            stress_constraint = d1_fea.beam_system.max_failure_crit * self.constraints_constants["stress_FOS"]
            constraints.append(stress_constraint)
            constraint_labels.append("stress")


        if active_constraints["energy"] or active_constraints["damage"] or active_constraints["natural frequency"] or active_constraints["frequency"]:
            def t_y(t):
                return self.total_force_y(t, drone_mass, drone_geometry.projected_surface_area)

            drone_power_module = DronePower(self.parameters["blade_num"], self.battery_system.milliwatt_hours,
                                            t_y, self.energy_calc_timesteps)
            drone_power_module.energy_consumption_calc()
            drone_power_module.throttle_ratio_trans_calc(16)

        freq_forcing_function = 0
        if active_constraints["natural frequency"] or active_constraints["frequency"]:
            propeller_RPM = self.parameters["propeller_max_RPM"] * max(drone_power_module.percent_throttle)
            freq_forcing_function = propeller_RPM * uc.rpm_to_rad_per_s * 2

        # Natural frequency constraint
        if active_constraints["natural frequency"]:
            d1_fea.solve_natural_frequencies()
            first_nf = d1_fea.beam_system.natural_frequencies[0]
            print("Propeller Freq",freq_forcing_function)
            print("NF:", d1_fea.beam_system.natural_frequencies)
            FOS_nf = first_nf / freq_forcing_function
            nf_constraint = self.constraints_constants["allowable_natural_frequency_FOS"] / FOS_nf

            constraints.append(nf_constraint)
            constraint_labels.append("natural frequency")


        if active_constraints["frequency"]:
            max_disp = max(d1_fea.beam_system.mag_displacements)
            d1_fea.beam_system.init_dynamic_forces(freq_forcing_function)
            d1_fea.beam_system.add_dynamic_harmonic_force(np.array([0,0,static_force,0,0,0]),"outer_node")
            d1_fea.solve_dynamic_harmonic(0,0)
            max_amp = max(d1_fea.beam_system.mag_displacements_amplitude)
            amp_disp_rat = max_amp / max_disp
            freq_constraint = amp_disp_rat / self.constraints_constants["allowable_amplitude_to_disp_ratio"]
            constraints.append(freq_constraint)
            constraint_labels.append("frequency")

        if active_constraints["energy"]:
            energy_used = drone_power_module.milliwatt_second_capacity - drone_power_module.energy_remaining[-1]

            energy_total = drone_power_module.milliwatt_second_capacity

            # Starts at 0 ends at 1
            energy_used_percent = energy_used / energy_total

            energy_constraint = energy_used_percent + self.constraints_constants["energy_remaining"]

            constraints.append(energy_constraint)
            constraint_labels.append("energy")

        if active_constraints["damage"]:
            surface_area = drone_geometry.projected_surface_area

            def f_t(t):
                return self.force_transient(t, drone_mass, surface_area)
            def f_o_t(t):
                return np.array([0, 0, self.force_oscillation_transient(t, drone_power_module, drone_mass, surface_area), 0, 0, 0])


            q = self.total_mission_duration / self.simulation_settings["transient_timesteps"]
            d1_fea.beam_system.initialize_dynamic_transient(q ,self.simulation_settings["transient_timesteps"], 0.25, 0.5)
            d1_fea.beam_system.add_dynamic_transient_force(f_t, "outer_node")
            d1_fea.beam_system.add_dynamic_transient_force(f_o_t, "outer_node")
            d1_fea.beam_system.solve_dynamic_transient()

            stress_curve = d1_fea.beam_system.beams[0].bending_stresses_0_ts
            detector = RF.FourPointDetector(recorder=RF.LoopValueRecorder())
            detector.process(stress_curve)
            collective = detector.recorder.collective
            cl = collective.load_collective
            damage_per_mission = self.woehler_curve_data_carbon_fiber.fatigue.damage(cl).sum()
            damage_constraint = damage_per_mission * self.constraints_constants["mission_count_damage"]

            constraints.append(damage_constraint)
            constraint_labels.append("damage")

        if detailed_output:
            # Nodes of intrest
            strut_node = d1_fea.beam_system.select_node('strut_node')
            strut_node_top = d1_fea.beam_system.select_node('strut_node_top')
            outer_node = d1_fea.beam_system.select_node('outer_node')
            center_node = d1_fea.beam_system.select_node('hub_node_0')

            # Elements of intrest
            hub_beam_0 = d1_fea.beam_system.select_element('hub_beam_0')
            strut_top_beam = d1_fea.beam_system.select_element("strut_element_top")
            outer_beam = d1_fea.beam_system.select_element("outer_beam")
            center_beam = d1_fea.beam_system.select_element("center_beam")

            strut_disp = d1_fea.beam_system.mag_displacements[strut_node]
            stress_strut = d1_fea.beam_system.beams[strut_top_beam].bending_stresses_0
            stress_hub = d1_fea.beam_system.beams[hub_beam_0].bending_stresses_0
            stress_outer = d1_fea.beam_system.beams[outer_beam].bending_stresses_0
            print("=== Optimizer Results ===")
            print(f"Wingspan: {wingspan} mm")
            print(f"Force: {-self.parameters["gravity"] * (self.parameters["payload"] / self.parameters["blade_num"])} N")
            print(f"Hub Beam Width: {(2/3) * hub_side_length} mm")
            print(f"Tip Displacement: {max_disp} mm")
            print(f"Strut Center Displacement: {strut_disp} mm")

            print(f"Stress Strut: {stress_strut} MPa")
            print(f"Stress Hub: {stress_hub} MPa")
            print(f"Stress Outer Beam: {stress_outer} MPa")

            print(f"Tip Displacement Amplitude: {max_amp} mm")
            print(f"Propeller Frequency {freq_forcing_function} rad/s")

        # When done LOG everything, then return constraints
        return np.array(constraints), constraint_labels

    def run_opt(self):
        dvm_list = np.array(list(self.design_variables_multipliers.values()))

        c_len = sum(list(self.active_constraints.values()))

        x_init = np.array(list(self.design_variables_initial_guess.values()))

        x_init_norm = x_init / dvm_list

        def mass_opt(x):
            d_v = {
                "arm_diameter": x[0] * dvm_list[0],
                "arm_thickness": x[1] * dvm_list[1],
                "strut_diameter": x[2] * dvm_list[2],
                "strut_thickness": x[3] * dvm_list[3],
                "strut_distance": x[4] * dvm_list[4],
                "hub_radius": x[5] * dvm_list[5],
                "hub_flange_thickness": x[6] * dvm_list[6],
                "hub_web_thickness": x[7] * dvm_list[7],
            }
            return self.mass(d_v) * self.simulation_settings["mass_multiplier"]
        def constraint_calculations_opt(x):

            d_v =  {
                "arm_diameter": x[0] * dvm_list[0],
                "arm_thickness": x[1] * dvm_list[1],
                "strut_diameter": x[2] * dvm_list[2],
                "strut_thickness": x[3] * dvm_list[3],
                "strut_distance": x[4] * dvm_list[4],
                "hub_radius": x[5] * dvm_list[5],
                "hub_flange_thickness": x[6] * dvm_list[6],
                "hub_web_thickness": x[7] * dvm_list[7],
            }
            res, labels = self.constraint_calculations(d_v, self.active_constraints)
            return res
        # Change to be an array of FOS

        # x_i' m(x') c(x') b' -> x_r' => x_r
        consts = (
            scipy.optimize.NonlinearConstraint(constraint_calculations_opt, np.zeros(c_len), np.ones(c_len))
        )

        boundaries = list(self.boundaries.values())
        boundaries_norm = []
        for i, bound in enumerate(boundaries):
            boundaries_norm.append((bound[0] / dvm_list[i], bound[1] / dvm_list[i]))
        results = scipy.optimize.minimize(mass_opt, x_init_norm, constraints=consts, bounds=boundaries_norm, method='SLSQP')

        return results

    def test_results(self, x, detailed_output=False):

        d_v = {
            "arm_diameter": x[0],
            "arm_thickness": x[1],
            "strut_diameter": x[2],
            "strut_thickness": x[3],
            "strut_distance": x[4],
            "hub_radius": x[5],
            "hub_flange_thickness": x[6],
            "hub_web_thickness": x[7],
        }
        return self.constraint_calculations(d_v, self.active_constraints_defaults,detailed_output=detailed_output)





