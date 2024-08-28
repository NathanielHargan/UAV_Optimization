#%%

#%%
import numpy as np
import matplotlib.pyplot as plt

class mission_profile:
    def __init__(self, name, initial_payload, acc_x, acc_y, max_vel_x, max_vel_y):
        # Name of the segment
        self.name = name

        # The ends of each segment are represented by values in the array. 0 for the start of the mission profile.
        self.t_values = np.array([0])
        self.x_coords = np.array([0])
        self.y_coords = np.array([0])

        self.acc_x = acc_x
        self.acc_y = acc_y
        self.max_vel_x = max_vel_x
        self.max_vel_y = max_vel_y

        self.payload = np.array([initial_payload])
        self.segment_names = ['origin']

        self.vel_x_values = np.array([0])
        self.vel_y_values = np.array([0])

        self.acc_x_values = np.array([]) # step function has 1 less value
        self.acc_y_values = np.array([])

    # Adds a new segment to the mission profile
    def add_segment_acc_y(self, name, x_end, y_end, y_vel_end):
        self.segment_names.append(name)
        self.segment_names.append(name)
        self.segment_names.append(name)

        self.payload = np.append(self.payload, self.payload[-1])

        y_start = self.y_coords[-1]
        y_v_start = self.vel_y_values[-1]

        # Segment 1
        t_1 = (self.max_vel_y - y_v_start) / self.acc_y  # time to reach max velocity y
        y_disp_1 = t_1 * y_v_start  + (self.acc_x * t_1 ** 2) / 2

        # Segment 3
        t_3 = (self.max_vel_y - y_vel_end) / self.acc_y
        y_disp_3 = t_3 * y_vel_end + (self.acc_x * t_3 ** 2) / 2

        # Segment 2
        y_disp_2 = y_end - y_start - y_disp_1 - y_disp_3
        t_2 = y_disp_2 / self.max_vel_y

        # Total time
        t_total = t_1 + t_2 + t_3

        # X change
        x_total_disp = x_end - self.vel_x_values[-1]
        x_velocity = x_total_disp/t_total

        x_disp_1 = x_velocity * t_1
        x_disp_2 = x_velocity * t_2
        x_disp_3 = x_velocity * t_3

        # append segment 1
        self.t_values = np.append(self.t_values, self.t_values[-1] + t_1)
        self.x_coords = np.append(self.x_coords, self.x_coords[-1] + x_disp_1)
        self.y_coords = np.append(self.y_coords, self.y_coords[-1] + y_disp_1)
        self.vel_x_values = np.append(self.vel_x_values, x_velocity)
        self.vel_y_values = np.append(self.vel_y_values, self.max_vel_y)
        self.acc_x_values = np.append(self.acc_x_values, 0)
        if y_end > y_start:
            self.acc_y_values = np.append(self.acc_y_values, self.acc_y)
        else:
            self.acc_y_values = np.append(self.acc_y_values, -self.acc_y)

        # append segment 2
        self.t_values = np.append(self.t_values, self.t_values[-1] + t_2)
        self.x_coords = np.append(self.x_coords, self.x_coords[-1] + x_disp_2)
        self.y_coords = np.append(self.y_coords, self.y_coords[-1] + y_disp_2)
        self.vel_x_values = np.append(self.vel_x_values, x_velocity)
        self.vel_y_values = np.append(self.vel_y_values, self.max_vel_y)
        self.acc_x_values = np.append(self.acc_x_values, 0)
        self.acc_y_values = np.append(self.acc_y_values, 0)

        # append segment 3
        self.t_values = np.append(self.t_values, self.t_values[-1] + t_3)
        self.x_coords = np.append(self.x_coords, self.x_coords[-1] + x_disp_3)
        self.y_coords = np.append(self.y_coords, self.y_coords[-1] + y_disp_3)
        self.vel_x_values = np.append(self.vel_x_values, x_velocity)
        self.vel_y_values = np.append(self.vel_y_values, y_vel_end)
        self.acc_x_values = np.append(self.acc_x_values, 0)
        if y_end > y_start:
            self.acc_y_values = np.append(self.acc_y_values, -self.acc_y)
        else:
            self.acc_y_values = np.append(self.acc_y_values, self.acc_y)

    # Adds a new segment to the mission profile
    def add_segment_acc_x(self, name, x_end, y_end, x_vel_end):

        self.segment_names.append(name)
        self.segment_names.append(name)
        self.segment_names.append(name)

        self.payload = np.append(self.payload, self.payload[-1])

        x_start = self.x_coords[-1]
        x_v_start = self.vel_x_values[-1]

        # Segment 1
        t_1 = (self.max_vel_x - x_v_start) / self.acc_x  # time to reach max velocity y
        x_disp_1 = t_1 * x_v_start + (self.acc_x * t_1 ** 2) / 2

        # Segment 3
        t_3 = (self.max_vel_y - x_vel_end) / self.acc_x
        x_disp_3 = t_3 * x_vel_end + (self.acc_x * t_3 ** 2) / 2

        # Segment 2
        x_disp_2 = x_end - x_start - x_disp_1 - x_disp_3
        t_2 = x_disp_2 / self.max_vel_x

        # Total time
        t_total = t_1 + t_2 + t_3

        # X change
        y_total_disp = y_end - self.vel_y_values[-1]
        y_velocity = y_total_disp/t_total

        y_disp_1 = y_velocity * t_1
        y_disp_2 = y_velocity * t_2
        y_disp_3 = y_velocity * t_3

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

        # append segment 2
        self.t_values = np.append(self.t_values, self.t_values[-1] + t_2)
        self.x_coords = np.append(self.x_coords, self.x_coords[-1] + x_disp_2)
        self.y_coords = np.append(self.y_coords, self.y_coords[-1] + y_disp_2)
        self.vel_x_values = np.append(self.vel_x_values, self.max_vel_x)
        self.vel_y_values = np.append(self.vel_y_values, y_velocity)
        self.acc_x_values = np.append(self.acc_x_values, 0)
        self.acc_y_values = np.append(self.acc_y_values, 0)

        # append segment 3
        self.t_values = np.append(self.t_values, self.t_values[-1] + t_3)
        self.x_coords = np.append(self.x_coords, self.x_coords[-1] + x_disp_3)
        self.y_coords = np.append(self.y_coords, self.y_coords[-1] + y_disp_3)
        self.vel_x_values = np.append(self.vel_x_values, x_vel_end)
        self.vel_y_values = np.append(self.vel_y_values, y_velocity)
        self.acc_y_values = np.append(self.acc_y_values, 0)
        if x_end > x_start:
            self.acc_x_values = np.append(self.acc_x_values, -self.acc_x)
        else:
            self.acc_x_values = np.append(self.acc_x_values, self.acc_x)

    def drop_payload(self, mass):
        self.payload[-1] = self.payload[-1] - mass

    # Calculates segment index and returns the position, velocity, and the name of the segment at a given time
    def time_solve(self, time):
        i = np.searchsorted(self.t_values, time, side='right') - 1  # Determines the segment time is in
        return self.time_output(time, i)

    # For a given value of time and segment index, it interpolates values for position, velocity, and the name of the segment
    def time_output(self, time, i):


        # Finds index of segment time is in
        # Side = 'right' means that if t == self.t_values[value], then it prioritizes the right one.

        t = (time - self.t_values[i])
        # Calculates time ratios
        delta_t = (self.t_values[i+1] - self.t_values[i])

        # Calculates position

        a_x = self.acc_x_values[i] + (t / delta_t) * (self.acc_x_values[i+1] - self.acc_x_values[i])
        a_y = self.acc_y_values[i] + (t / delta_t) * (self.acc_y_values[i+1] - self.acc_y_values[i])

        v_x = self.vel_x_values[i] + t * self.acc_x_values[i] + ((t ** 2) / (2 * delta_t)) * (self.acc_x_values[i+1] - self.acc_x_values[i])
        v_y = self.vel_y_values[i] + t * self.acc_y_values[i] + ((t ** 2) / (2 * delta_t)) * (self.acc_y_values[i+1] - self.acc_y_values[i])

        x = self.x_coords[i] + t * self.vel_x_values[i] + (t ** 2 / 2) * self.acc_x_values[i] + ((t ** 3) / (6 * delta_t)) * (self.acc_x_values[i+1] - self.acc_x_values[i])
        y = self.y_coords[i] + t * self.vel_y_values[i] + (t ** 2 / 2) * self.acc_y_values[i] + ((t ** 3) / (6 * delta_t)) * (self.acc_y_values[i+1] - self.acc_y_values[i])

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
    print("Y acc:")
    print(np.round(m1.acc_y_values, 2))

    list_y = np.empty(0)
    list_t = np.empty(0)
    list_n = []
    for t in range(200):
        x, y, v_x, v_y, a_x, a_y, name = m1.time_solve(t)
        print(name + " - y: " + str(y))

if __name__ == '__main__':
    test_mission_profile()
