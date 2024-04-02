import numpy as np
class mission_profile:
    def __init__(self, name, initial_payload):
        # Name of the segment
        self.name = name

        # The ends of each segment are represented by values in the array. 0 for the start of the mission profile.
        self.t_values = np.array([0])
        self.x_coords = np.array([0])
        self.y_coords = np.array([0])
        self.payload = np.array([initial_payload])

        # One velocity for each segment
        self.v_x_values = np.array([])
        self.v_y_values = np.array([])

        self.segment_names = []

        self.thrust = np.array([])

    # Adds a new segment to the mission profile
    def add_segment(self,name,t,x_end=None,y_end=None):

        # Adds segment name to array
        self.segment_names.append(name)

        self.payload = np.append(self.payload, self.payload[-1])

        # If X and Y coordinates are not input it assumes loitering.
        if x_end is None:
            self.x_coords = np.append(self.x_coords, self.x_coords[-1]) # Repeat the last X coordinate
            self.y_coords = np.append(self.y_coords, self.y_coords[-1]) # Repeat the last Y coordinate
            self.v_x_values = np.append(self.v_x_values, 0) # Adds 0 to velocity
            self.v_y_values = np.append(self.v_y_values, 0) # Adds 0 to velocity

        else:
            self.x_coords = np.append(self.x_coords, x_end) # New X is added
            self.y_coords = np.append(self.y_coords, y_end) # New Y is added

            vx = (self.x_coords[-1] - self.x_coords[-2])/t # Calculates the velocity in X
            vy = (self.y_coords[-1] - self.y_coords[-2])/t # Calculates the velocity in Y

            self.v_x_values = np.append(self.v_x_values, vx) # Appends the velocity in X
            self.v_y_values = np.append(self.v_y_values, vy) # Appends the velocity in Y

        self.t_values = np.append(self.t_values,t + self.t_values[-1]) # Adds time to last time recorded.

    def drop_payload(self, mass):
        self.payload[-1] = self.payload[-1] - mass

    # Calculates segment index and returns the position, velocity, and the name of the segment at a given time
    def time_solve(self,time):
        i = np.searchsorted(self.t_values, time, side='right') - 1 # Determines the segment time is in
        return self.time_output(time,i)

    # For a given value of time and segment index, it interpolates values for position, velocity, and the name of the segment
    def time_output(self,time,i):
        # Finds index of segment time is in
        # Side = 'right' means that if t == self.t_values[value], then it prioritizes the right one.

        # Calculates time ratios
        n1 = (self.t_values[i+1] - time) / (self.t_values[i+1] - self.t_values[i])
        n2 = (time - self.t_values[i]) / (self.t_values[i+1] - self.t_values[i])

        # Calculates position
        x = n1 * self.x_coords[i] + n2 * self.x_coords[i+1]
        y = n1 * self.y_coords[i] + n2 * self.y_coords[i+1]

        v_x = self.v_x_values[i]
        v_y = self.v_x_values[i]

        yield x # Outputs position X
        yield y # Outputs position Y
        yield v_x # Outputs velocity X
        yield v_y # Outputs velocity Y
        yield self.segment_names[i] # Outputs segment name
        return # Ends the function


#%%
