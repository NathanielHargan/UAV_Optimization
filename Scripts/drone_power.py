from Scripts.material_properties import Gravity
import numpy as np


class BatterySystem:
    def __init__(self, name, volts_battery, amp_hours_battery, battery_series_num):
        self.name = name
        self.volts_battery = volts_battery
        self.milliamp_hours_battery = amp_hours_battery
        self.battery_series_num = battery_series_num
        self.volts_total = volts_battery * battery_series_num
        self.milliwatt_hours = amp_hours_battery * self.volts_total


class DronePower:
    def __init__(self, blade_num, mwh, drone_forces):
        self.blade_num = blade_num
        self.milliwatt_hour_capacity = mwh
        self.milliwatt_second_capacity = mwh * (3600)
        self.drone_forces = drone_forces

    def thrust_to_power(self, thrust):
        # The plot given looks like a parabola with points intersecting at:
        # (thrust,power): (10kg,1kW)/8 , (50kg, 6kW)/8 , (80kg, 13kW)/8
        thrust_mass_per_blade_kg = 1000 * thrust / (Gravity * self.blade_num)  # N => N / (mm/s^2) * (1000 mm/m) = kg
        power_kilowatts_per_blade = (0.01238 * thrust_mass_per_blade_kg ** 2 + 0.03214 * thrust_mass_per_blade_kg + 0.06548)
        power_milliwatts_per_blade = power_kilowatts_per_blade * 10 ** 6

        return power_milliwatts_per_blade * self.blade_num

    def power_at_time(self, t):
        thrust_y, drag = self.drone_forces.time_solve(t)
        power = self.thrust_to_power(thrust_y)

        return power

    def energy_consumption(self, steps, segments):
        # Segment ends are a list of the segment deviations.
        segment_count = len(segments)-1
        time = np.empty([])

        # Setting up time segments
        for i in range(segment_count):
            t_start = segments[i]
            t_end = segments[i+1]
            time = np.append(time, np.arange(t_start, t_end, steps))

        power_consumed = np.array([self.milliwatt_second_capacity])  # Initial Power
        for i in range(len(time)-1):
            # Runge Kutta implementation
            t0 = time[i]
            t1 = time[i+1]
            h = t1 - t0

            # negative because power used results in a loss of energy
            k1 = -self.power_at_time(t0)
            k2 = -self.power_at_time(t0 + h/2)
            k3 = -self.power_at_time(t0 + h/2)
            k4 = -self.power_at_time(t0 + h)

            next_power_consumed = power_consumed[-1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

            power_consumed = np.append(power_consumed, next_power_consumed)

        yield time  # Array from 0 to t_final
        yield power_consumed
        return

if __name__ == '__main__':
    from Scripts.drone_mass import DroneMass
    from Scripts.drone_battery import Battery
    from Scripts.drone_geometry import DroneGeometry
    from Scripts.beams import BeamType
    from Scripts.mission_profile import MissionProfile

    # Mass Analysis Setup
    battery_1 = Battery(4.292e-3, np.array([0, 86.75, -50]), np.array([260.5, 123.5, 63.5]))
    battery_2 = Battery(4.292e-3, np.array([0, -86.75, -50]), np.array([260.5, 123.5, 63.5]))
    arm_beam = BeamType("annulus", [29, 25], "Aluminum7075-T6", "arm_beam")
    strut_beam = BeamType("annulus", [14, 12], "Aluminum7075-T6", "strut_beam")
    d1 = DroneGeometry("drone_test", 1000, 660, 8, arm_beam, strut_beam, [battery_1, battery_2])
    d1_mass = DroneMass(d1)

    # Mass Analysis results
    total_mass = d1_mass.total_mass
    surface_area = d1.projected_surface_area

    # Mission Profile generation
    m1 = MissionProfile("m1", 1, 1, 0)
    m1.add_segment("takeoff",10,50,5,5,1)
    m1.add_segment("climb",30,550,1000,30,1)
    m1.add_segment("cruise",100,1000,1100,30,1)
    m1.add_segment("descent",30,1100,500,0,0)
    m1.add_segment("loiter",30,1100,500,0,0)
    m1.add_segment("land",20,1200,0,0,0)

    # Drone Force Analysis
    d1_forces = DroneForces(total_mass, surface_area, 1, m1)

    # Drone
    battery_sys = BatterySystem(
    "Tattu 40000mAh 6S 10C 22.8V High Voltage UAV Lipo Battery Pack with AS150+AS150",
    22.8,
    40,
    2)
    drone_power = DronePower(8, battery_sys.milliwatt_hours, d1_forces)


#%%
