#%%

#%%
import numpy as np
import scipy
import matplotlib.pyplot as plt


class MissionProfileCubicSpline:
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
        if time < 0 or time > self.t_values[-1]:
            return np.array([0, 0, 0, 0]), np.array([0, 0, 0, 0]), "NA"
        else:
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

def kinematic_matrix(dt):
    # [px0, py0, vx0, vy0, ax0, ay0, px1, py1, vx1, vy1, ax1, ay1]
    k_m = np.array([
        [1, 0, dt, 0, (1/3)*dt**2, 0, -1, 0, 0, 0, (1/6)*dt**2, 0],    # x
        [0, 1, 0, dt, 0, (1/3)*dt**2, 0, -1, 0, 0, 0, (1/6)*dt**2],    # y
        [0, 0, 1, 0, (1/2)*dt, 0, 0, 0, -1, 0, (1/2)*dt, 0],           # x'
        [0, 0, 0, 1, 0, (1/2)*dt, 0, 0, 0, -1, 0, (1/2)*dt],           # y'
    ])

    return k_m

class MissionProfileSegment:
    def __init__(self, dt, name='node'):
        self.name = name
        self.dt = dt
        self.local_kinematic_matrix = kinematic_matrix(dt)

class MissionProfileLinearAcc:
    def __init__(self, name='profile'):
        self.name = name
        self.nodes = np.array([0])
        self.node_names = ['origin']

        self.segment_names = []
        self.segments = []

        self.constraints_indexes = np.array([])
        self.constraints_values = np.array([])
        self.duration = 0

    def add_segment(self, dt, node_name='node', seg_name='segment'):

        if node_name == 'node':
            node_name = 'node'+str(len(self.nodes))

        if seg_name == 'segment':
            seg_name = 'segment'+str(len(self.segments))

        self.segments.append(MissionProfileSegment(dt, seg_name))
        self.segment_names.append(seg_name)


        self.nodes = np.append(self.nodes, self.nodes[-1] + dt)
        self.node_names.append(node_name)



    def constraint_index(self, ref):
        if type(ref) is str:
            return {"x": 0, "y": 1, "vx": 2, "vy": 3, "ax": 4, "ay": 5}[ref]
        elif type(ref) is int:
            return ref
        else:
            print("error")
            return "error"

    def add_constraint(self, value, ref, node):
        cons_index_local = self.constraint_index(ref)
        cons_index_global = int(6 * node + cons_index_local)
        self.constraints_indexes = np.append(self.constraints_indexes, cons_index_global)
        self.constraints_values = np.append(self.constraints_values, value)

    def solve_kinematics(self):
        node_num = len(self.segments) + 1
        global_kinematic_matrix = np.zeros([4 * node_num, 6 * node_num])
        # Assemble the matrix
        for i, segment in enumerate(self.segments):
            global_kinematic_matrix[4*i:4*i+4,6*i:6*i+12] += segment.local_kinematic_matrix

        u_ind = self.constraints_indexes.astype(int)
        p_ind = np.delete(np.arange(0, node_num * 6), u_ind)
        km_p = np.take(global_kinematic_matrix[:], self.constraints_indexes.astype(int), axis=1)
        km_u = np.delete(global_kinematic_matrix[:], self.constraints_indexes.astype(int), axis=1)

        b = -km_p @ self.constraints_values
        # for segment in self.segments:
        solution = np.linalg.lstsq(km_u, b, rcond=None)[0]

        solution_vector = np.zeros(node_num*6)

        for i, u_i in enumerate(u_ind):
            solution_vector[u_i] = self.constraints_values[i]

        for i, p_i in enumerate(p_ind):
            solution_vector[p_i] = solution[i]

        self.solution_vector = solution_vector

        self.px = solution_vector[0::6]
        self.py = solution_vector[1::6]
        self.vx = solution_vector[2::6]
        self.vy = solution_vector[3::6]
        self.ax = solution_vector[4::6]
        self.ay = solution_vector[5::6]

        self.duration = self.nodes[-1]

    def time_solve(self, t):
        # find index
        seg_i = np.searchsorted(self.nodes, t, side='right') - 1
        dt = self.segments[seg_i].dt
        tau = t - self.nodes[seg_i]
        ax = self.ax[seg_i] + (tau/dt) * (self.ax[seg_i+1] - self.ax[seg_i])
        vx = self.vx[seg_i] + tau * self.ax[seg_i] + (tau**2 / (2*dt)) * (self.ax[seg_i+1] - self.ax[seg_i])
        px = self.px[seg_i] + tau * self.vx[seg_i] + (tau**2 / 2) * self.ax[seg_i] + (tau**3 / (6*dt)) * (self.ax[seg_i+1] - self.ax[seg_i])

        ay = self.ay[seg_i] + (tau/dt) * (self.ay[seg_i+1] - self.ay[seg_i])
        vy = self.vy[seg_i] + tau * self.ay[seg_i] + (tau**2 / (2*dt)) * (self.ay[seg_i+1] - self.ay[seg_i])
        py = self.py[seg_i] + tau * self.vy[seg_i] + (tau**2 / 2) * self.ay[seg_i] + (tau**3 / (6*dt)) * (self.ay[seg_i+1] - self.ay[seg_i])

        return np.array([px,py]), np.array([vx,vy]), np.array([ax,ay])

def mission_profile():

    m1 = MissionProfileLinearAcc()

    m1.add_segment(15, "takeoff")
    m1.add_segment(15, "takeoff2")
    m1.add_constraint(0, "x", 0)
    m1.add_constraint(0, "y", 0)
    m1.add_constraint(0, "vx", 0)
    m1.add_constraint(0, "vy", 0)
    m1.add_constraint(0, "ax", 0)
    m1.add_constraint(0, "ay", 0)
    m1.add_constraint(10, "x", 2)
    m1.add_constraint(20, "y", 2)
    m1.add_constraint(0, "vx", 2)
    m1.add_constraint(0, "vy", 2)
    m1.add_constraint(0, "ax", 2)
    m1.add_constraint(0, "ay", 2)
    m1.solve_kinematics()
    print(m1.solution_vector)




if __name__ == '__main__':
    mission_profile()

#%%
