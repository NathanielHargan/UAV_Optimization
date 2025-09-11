import scipy
import numpy as np
import Scripts.unit_conversions as uc

class BatterySystem:
    def __init__(self, name, volts_battery, milliamp_hours_battery, battery_series_num):
        self.name = name
        self.volts_battery = volts_battery
        self.milliamp_hours_battery = milliamp_hours_battery
        self.battery_series_num = battery_series_num
        self.volts_total = volts_battery * battery_series_num
        self.milliwatt_hours = milliamp_hours_battery * self.volts_total


class DronePower:
    def __init__(self, blade_num, mwh, total_thrust_y, timesteps):
        self.blade_num = blade_num
        self.milliwatt_hour_capacity = mwh
        self.milliwatt_second_capacity = mwh * (uc.hrs_to_s)
        self.drone_forces = total_thrust_y
        self.energy_remaining = np.array([])
        self.frequency = np.array([])
        self.ratio_throttle = np.array([])
        self.error = np.array([])
        self.timesteps = timesteps

    def thrust_to_power(self, thrust):
        # (thrust,power): (10kg,1kW)/8 , (50kg, 6kW)/8 , (80kg, 13kW)/8

        thrust_mass_per_blade_kg = uc.N_to_kgf * thrust / self.blade_num  # N => (N / (mm/s^2)) * (1000 mm/m) = kg
        power_kilowatts_per_blade = (0.01238 * thrust_mass_per_blade_kg ** 2 + 0.03214 * thrust_mass_per_blade_kg + 0.06548)
        power_milliwatts_per_blade = power_kilowatts_per_blade * uc.kW_to_mW

        return power_milliwatts_per_blade * self.blade_num

    def power_at_time(self, t):
        thrust_y = self.drone_forces(t)
        power = self.thrust_to_power(thrust_y)

        return power

    def neg_power_at_time(self, t):
        return -self.power_at_time(t)

    '''
    def energy_consumption(self, t_end):
        integration = scipy.integrate.solve_ivp(self.neg_power_at_time, [0.1, t_end], [self.milliwatt_second_capacity], method = 'DOP853')
        yield integration.t
        yield integration.y
        return
    '''
    '''
    def energy_consumption_calc(self, segments, steps_mult, steps_init, max_error):
        # Segment ends are a list of the segment deviations.
        segment_count = len(segments)-1
        steps = steps_init * steps_mult

        power_consumed_old = np.ones(steps_init * segment_count) * self.milliwatt_second_capacity
        power_consumed = np.array([])
        error = np.ones(steps_init * segment_count) * self.milliwatt_second_capacity

        while max(error) > max_error:
            time = np.array([])
            
            # Setting up time segments
            for i in range(segment_count):
                t_start = segments[i]
                t_end = segments[i+1]
                time = np.append(time, np.linspace(t_start, t_end, steps, endpoint = False))

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

            error = np.abs((power_consumed[::steps_mult] / power_consumed_old) - 1)
            steps = steps * steps_mult
            power_consumed_old = power_consumed

        print("Steps taken: ", steps)
        print("Max error: ", max(error))
        self.error = error
        self.timesteps = time  # Array from 0 to t_final
        self.energy_consumption = power_consumed
    '''
    def energy_consumption_calc(self):
        power = np.zeros(len(self.timesteps))
        for i, t in enumerate(self.timesteps):
            power[i] = self.neg_power_at_time(float(t))
        self.power = power
        self.energy_remaining = self.milliwatt_second_capacity + scipy.integrate.cumulative_trapezoid(power, self.timesteps, initial = 0)

    def throttle_ratio_trans_calc(self, max_power_full_battery_kw):
        # assuming max throttle is proportional to the percent energy left in the battery.
        self.percent_throttle = np.zeros(len(self.timesteps))
        max_throttle_power = (max_power_full_battery_kw * uc.kW_to_mW)
        for i, t in enumerate(self.timesteps):
            # percent_energy = self.energy_consumption[i] / self.milliwatt_second_capacity # % energy in the battery
            # max_throttle_power = (max_power_full_battery_kw * uc.kW_to_mW) * percent_energy
            self.percent_throttle[i] = self.power_at_time(t) / (max_throttle_power)

    def throttle_ratio_calc(self, max_power_full_battery_kw, t):
        max_throttle_power = (max_power_full_battery_kw * uc.kW_to_mW)
        return self.power_at_time(t) / (max_throttle_power)

    def freq_trans_calc(self, rpm_max):
        frequency = np.zeros(len(self.timesteps))
        for i, t in enumerate(self.timesteps):
            frequency[i] = uc.Hz_to_rad_per_s * rpm_max * self.percent_throttle[i] / 30
        self.frequency = frequency

    def freq_calc(self, rpm_max, max_power_full_battery_kw, t):
        return uc.Hz_to_rad_per_s * rpm_max * self.throttle_ratio_calc(max_power_full_battery_kw, t) / 30

#%%
