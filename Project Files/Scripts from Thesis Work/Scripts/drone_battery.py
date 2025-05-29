import numpy as np


class Battery:
    def __init__(self, mass, position, dimensions):
        self.mass = mass
        self.position = position # position of the center of the battery
        self.dimensions = dimensions

        l = dimensions[0]
        h = dimensions[1]
        w = dimensions[2]

        self.local_moment_of_inertia = np.array([[(1/12) * mass * (h**2 + w**2), 0, 0],
                                                 [0, (1/12) * mass * (l**2 + w**2), 0],
                                                 [0, 0, (1/12) * mass * (l**2 + h**2)]])

        self.transformation = np.array([[1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]])
