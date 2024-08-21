#%%

#%%
import numpy as np
import matplotlib.pyplot as plt

class mission_profile:
    def __init__(self, name, initial_payload):
        # Name of the segment
        self.name = name

        # The ends of each segment are represented by values in the array. 0 for the start of the mission profile.
        self.t_values = np.array([0])
        self.x_coords = np.array([0])
        self.y_coords = np.array([0])
        self.payload = np.array([initial_payload])
        self.segment_names = []

        self.vel_x_values = np.array([0])
        self.vel_y_values = np.array([0])

        self.acc_x_values = np.array([0])
        self.acc_y_values = np.array([0])

    # Adds a new segment to the mission profile
    def add_segment(self, name, t, x_end=None, y_end=None, x_vel_end=None, y_vel_end=None, x_a_end=None, y_a_end=None):

        # Adds segment name to array
        self.segment_names.append(name)
        self.payload = np.append(self.payload, self.payload[-1])
        self.t_values = np.append(self.t_values, self.t_values[-1] + t)

        # If X and Y coordinates are not input it assumes loitering.
        if x_end is None:
            self.x_coords = np.append(self.x_coords, self.x_coords[-1]) # Repeat the last X coordinate
            self.y_coords = np.append(self.y_coords, self.y_coords[-1]) # Repeat the last Y coordinate
        else:
            self.x_coords = np.append(self.x_coords, x_end)
            self.y_coords = np.append(self.y_coords, y_end)

        if x_a_end is None:
            self.acc_x_values = np.append(self.acc_x_values,0)
            self.acc_y_values = np.append(self.acc_y_values,0)
        else:
            self.acc_x_values = np.append(self.acc_x_values, x_a_end)
            self.acc_y_values = np.append(self.acc_y_values, y_a_end)

        if x_vel_end is None:
            self.vel_x_values = np.append(self.vel_x_values,0)
            self.vel_y_values = np.append(self.vel_y_values,0)
        else:
            self.vel_x_values = np.append(self.vel_x_values, x_a_end)
            self.vel_y_values = np.append(self.vel_y_values, y_a_end)

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
        print(t)
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

    m1 = mission_profile("m1", 50)
    m1.add_segment("takeoff", 10, 50, 5, 1, 10, 0, 9.81*0.2)
    m1.add_segment("climb", 30, 550, 1000, 1, 1, 0, 0)
    m1.add_segment("cruise", 100, 1000, 11001, 1, 0, 0, 0)
    m1.add_segment("descent", 30, 1100, 5001, 1, 0, 0, 0)
    m1.add_segment("loiter", 30)
    m1.add_segment("land", 20, 0, 0)

    list_y = np.empty(0)
    list_t = np.empty(0)
    list_n = []
    for t in range(200):
        x, y, v_x, v_y, a_x, a_y, name = m1.time_solve(t)

        list_t = np.append(list_t, t)
        list_y = np.append(list_y, y)
        list_n.append(name)
    print(list_y)
    plt.plot(list_n, list_y)

if __name__ == '__main__':
    test_mission_profile()
