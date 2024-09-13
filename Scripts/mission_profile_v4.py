#%%

#%%
import numpy as np
import matplotlib.pyplot as plt

def shape_function_timoshenko_a_inv(L, g):
    a_inv = (1 / (L ** 2 + 12 * g)) * np.array([
        [L ** 2 + 12 * g, 0, 0, 0],
        [-12 / g, L ** 2 + 6 * g, 12 * g / L, -6 * g],
        [-3, -(2 * L ** 2 + 6 * g)/L, 3, -(L ** 2 - 6 * g)/L],
        [2/L, 1, -2/L, 1]
    ])
    return a_inv


class mission_profile:
    def __init__(self, name, initial_payload, initial_a_x, initial_a_y):
        # Name of the segment
        self.name = name

        # The ends of each segment are represented by values in the array. 0 for the start of the mission profile.
        self.t_values = np.array([0])
        self.x_coords = np.array([0])
        self.y_coords = np.array([0])

        self.payload = np.array([initial_payload])
        self.segment_names = ['origin']

        self.acc_x_values = np.array([initial_a_x]) # step function has 1 less value
        self.acc_y_values = np.array([initial_a_y])

    # Adds a new segment to the mission profile
    def add_segment(self, name, x_end, y_end, x_curve, y_curve):
        self.segment_names.append(name)
        self.payload = np.append(self.payload, self.payload[-1])


        # append segment 1
        self.t_values = np.append(self.t_values, self.t_values[-1] + t_1)
        self.x_coords = np.append(self.x_coords, self.x_coords[-1] + x_disp_1)
        self.y_coords = np.append(self.y_coords, self.y_coords[-1] + y_disp_1)
        self.vel_x_values = np.append(self.vel_x_values, self.max_vel_x)
        self.vel_y_values = np.append(self.vel_y_values, y_velocity)
        self.acc_y_values = np.append(self.acc_y_values, 0)
        if x_end > x_start:
            self.acc_x_values = np.append(self.acc_x_values, self.acc_x)
        else:
            self.acc_x_values = np.append(self.acc_x_values, -self.acc_x)



    def drop_payload(self, mass):
        self.payload[-1] = self.payload[-1] - mass

    # Calculates segment index and returns the position, velocity, and the name of the segment at a given time
    def time_solve(self, time):
        i = np.searchsorted(self.t_values, time, side='right') - 1  # Determines the segment time is in
        return self.time_output(time, i)

    # For a given value of time and segment index, it interpolates values for position, velocity, and the name of the segment
    def time_output(self, time, i):



        vector = np.array([self.x_coords[i],
                           self.acc_x_values[i],
                           self.y_coords[i],
                           self.acc_y_values[i],
                           self.x_coords[i+1],
                           self.acc_x_values[i+1],
                           self.y_coords[i+1],
                           self.acc_y_values[i+1]])

        timestep_transform = np.array([
            [0,0,0,0,0,0],
            [0,0,0,0,0,0],
            [0,0,0,0,0,0]])

        # Finds index of segment time is in
        # Side = 'right' means that if t == self.t_values[value], then it prioritizes the right one.

        t = (time - self.t_values[i])
        # Calculates time ratios
        delta_t = (self.t_values[i+1] - self.t_values[i])

        # Calculates position

        a_x = self.acc_x_values[i] + (t / delta_t) * (self.acc_x_values[i] - self.acc_x_values[i])
        a_y = self.acc_y_values[i] + (t / delta_t) * (self.acc_y_values[i] - self.acc_y_values[i])

        v_x = self.vel_x_values[i] + t * self.acc_x_values[i] + ((t ** 2) / (2 * delta_t)) * (self.acc_x_values[i] - self.acc_x_values[i])
        v_y = self.vel_y_values[i] + t * self.acc_y_values[i] + ((t ** 2) / (2 * delta_t)) * (self.acc_y_values[i] - self.acc_y_values[i])

        x = self.x_coords[i] + t * self.vel_x_values[i] + (t ** 2 / 2) * self.acc_x_values[i] + ((t ** 3) / (6 * delta_t)) * (self.acc_x_values[i] - self.acc_x_values[i])
        y = self.y_coords[i] + t * self.vel_y_values[i] + (t ** 2 / 2) * self.acc_y_values[i] + ((t ** 3) / (6 * delta_t)) * (self.acc_y_values[i] - self.acc_y_values[i])

        yield x # Outputs position X
        yield y # Outputs position Y
        yield v_x # Outputs velocity X
        yield v_y # Outputs velocity Y
        yield a_x
        yield a_y
        yield self.segment_names[i] # Outputs segment name
        return # Ends the function

def test_mission_profile():

    m1 = mission_profile("m1",
                         50,
                         4,
                         0.2*9.8, # acc_y
                         10, 10)

    m1.add_segment_acc_x("takeoff", 60, 10, 5)
    m1.add_segment_acc_y("climb", 60, 1000, 2)

    print("\nnames:")
    print(m1.segment_names)
    print("time:")
    print(np.round(m1.t_values,2))
    print("X:")
    print(np.round(m1.x_coords,2))
    print("Y:")
    print(np.round(m1.y_coords, 2))
    print("X vel:")
    print(np.round(m1.vel_x_values, 2))
    print("Y vel:")
    print(np.round(m1.vel_y_values, 2))
    print("X acc:")
    print(np.round(m1.acc_x_values, 2))
    print("Y acc:")
    print(np.round(m1.acc_y_values, 2))



if __name__ == '__main__':
    test_mission_profile()

#%%

        if len(self.vel_x_values) == 1:
            self.vel_x_values = np.array([x_velocity])
#%%
