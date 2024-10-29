#%%

#%%
import numpy as np
import matplotlib.pyplot as plt


class MissionProfile:
    def __init__(self, name, initial_payload, initial_v_x, initial_v_y):
        # Name of the segment
        self.name = name

        # The ends of each segment are represented by values in the array. 0 for the start of the mission profile.
        self.t_values = np.array([0])
        self.x_coords = np.array([0])
        self.y_coords = np.array([0])

        self.payload = np.array([initial_payload])
        self.segment_names = ['origin']

        self.vel_x_values = np.array([initial_v_x]) # step function has 1 less value
        self.vel_y_values = np.array([initial_v_y])

    # Adds a new segment to the mission profile
    def add_segment(self, name, time, x, y, v_x, v_y):
        self.segment_names.append(name)
        self.payload = np.append(self.payload, self.payload[-1])

        self.t_values = np.append(self.t_values, self.t_values[-1] + time)
        self.x_coords = np.append(self.x_coords, x)
        self.y_coords = np.append(self.y_coords, y)
        self.vel_y_values = np.append(self.vel_y_values, v_y)
        self.vel_x_values = np.append(self.vel_x_values, v_x)

    def drop_payload(self, mass):
        self.payload[-1] = self.payload[-1] - mass

    # Calculates segment index and returns the position, velocity, and the name of the segment at a given time
    def time_solve(self, time):
        i = np.searchsorted(self.t_values, time, side='right') - 1  # Determines the segment time is in
        return self.time_output(time, i)


    # For a given value of time and segment index, it interpolates values for position, velocity, and the name of the segment
    def time_output(self, time, i):
        dt = (self.t_values[i+1] - self.t_values[i])
        ti = (time - self.t_values[i]) / dt

        vector_y = np.array([self.y_coords[i],
                             self.vel_y_values[i] * dt,
                             self.y_coords[i+1],
                             self.vel_y_values[i+1] * dt])

        vector_x = np.array([self.x_coords[i],
                             self.vel_x_values[i] * dt,
                             self.x_coords[i+1],
                             self.vel_x_values[i+1] * dt])

        # https://en.wikipedia.org/wiki/Cubic_Hermite_spline
        time_matrix = np.array([[2*ti**3 - 3*ti**2 + 1, ti**3 - 2*ti**2 + ti, -2*ti**3 + 3*ti**2, ti**3 - ti**2],
                                [(6*ti**2 - 6*ti) / dt, (3*ti**2 - 4*ti + 1) / dt, (-6*ti**2 + 6*ti) / dt, (3*ti**2 - 2*ti) / dt],
                                [(12*ti - 6) / dt ** 2, (6*ti - 4)/ dt ** 2, (-12*ti + 6) / dt ** 2, (6*ti - 2) / dt ** 2],
                                [12 / dt ** 3, 6 / dt ** 3, -12 / dt ** 3, 6 / dt ** 3]])

        x_result = time_matrix @ np.transpose(vector_x)
        y_result = time_matrix @ np.transpose(vector_y)

        # Finds index of segment time is in
        # Side = 'right' means that if t == self.t_values[value], then it prioritizes the right one.

        yield x_result # Outputs position X
        yield y_result # Outputs position Y
        yield self.segment_names[i] # Outputs segment name
        return # Ends the function


def test_mission_profile():

    m1 = MissionProfile("m1",
                         50,
                         1,
                         2)

    m1.add_segment("takeoff", 60, 10, 5)
    m1.add_segment("climb", 60, 1000, 2)
    m1.add_segment("climb 2", 60, 1100, 1)

    print("\nnames:")
    print(m1.segment_names)
    print("time:")
    print(np.round(m1.t_values,2))
    print("X:")
    print(np.round(m1.x_coords,2))
    print("Y:")



if __name__ == '__main__':
    test_mission_profile()

#%%
#%%
