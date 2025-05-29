import numpy as np

def mass_moment_annulus(cross_section_parameters, L, density):
    r_2 = cross_section_parameters[0]/2  # outer
    r_1 = cross_section_parameters[1]/2  # inner
    m = np.pi * (r_2 ** 2 - r_1 ** 2) * L * density
    I = np.array([
        [0.5 * m * (r_2 ** 2 + r_2 ** 2), 0, 0],
        [0, (1/12) * m * (3 * (r_2 ** 2 + r_2 ** 2) ** 2 + L ** 2), 0],
        [0, 0, (1/12) * m * (3 * (r_2 ** 2 + r_2 ** 2) ** 2 + L ** 2)]
    ])
    return I


def mass_moment_circle(cross_section_parameters, L, density):
    r = cross_section_parameters[0]/2
    m = np.pi * r ** 2 * L * density
    I = np.array([
        [0.5 * m * r ** 2, 0, 0],
        [0, (1/12) * m * (3 * r ** 2 + L ** 2), 0],
        [0, 0, (1/12) * m * (3 * r ** 2 + L ** 2)]
    ])
    return I

def mass_moment_rectangle(cross_section_parameters, L, density):
    h = cross_section_parameters[0]
    b = cross_section_parameters[1]
    m = h * b * L * density
    I = np.array([
        [(1/12) * m * (h ** 2 + b ** 2), 0, 0],
        [0, (1/12) * m * (h ** 2 + b ** 2), 0],
        [0, 0, (1/12) * m * (h ** 2 + b ** 2)]
    ])
    return I

def mass_moment_hexagon(cross_section_parameters, L, density):
    h = cross_section_parameters[0]
    b = cross_section_parameters[1]
    m = h * b * L * density
    I = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ])
    return I

