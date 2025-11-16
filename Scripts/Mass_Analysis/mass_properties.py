import numpy as np
import math
def mass_moment_annulus(cross_section_parameters, L, density):
    r_2 = cross_section_parameters[0]/2  # outer
    r_1 = cross_section_parameters[1]/2  # inner
    m = np.pi * (r_2 ** 2 - r_1 ** 2) * L * density
    I = np.array([
        [0.5 * m * (r_2 ** 2 + r_1 ** 2), 0, 0],
        [0, (1/12) * m * (3 * (r_2 ** 2 + r_1 ** 2) ** 2 + L ** 2), 0],
        [0, 0, (1/12) * m * (3 * (r_2 ** 2 + r_1 ** 2) ** 2 + L ** 2)]
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
        [0, (1/12) * m * (h ** 2 + L ** 2), 0],
        [0, 0, (1/12) * m * (L ** 2 + b ** 2)]
    ])
    return I

def mass_moment_hexagon(cross_section_parameters, L, density):
    r = cross_section_parameters[0]
    m = (3/8) * math.sqrt(3) * (r ** 2) * L * density
    I = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ])
    return I

def mass_moment_i_beam(cross_section_parameters, L, density):
    b = cross_section_parameters[0]
    d = cross_section_parameters[1]
    s = cross_section_parameters[2]
    t = cross_section_parameters[3]

    h = d - 2 * s
    cz = b / 2
    cy = d / 2
    # Area
    area = b * d - h * (b - t)
    m = area * L * density
    m_t = s * b * L * density# mass of top
    m_w = h * t * L * density# mass of web


    # Centroid I beam.
    I_top = np.array([
        [(1 / 12) * m_t * (s ** 2 + b ** 2), 0, 0],
        [0, (1 / 12) * m_t * (L ** 2 + b ** 2), 0],
        [0, 0, (1 / 12) * m_t * (s ** 2 + L ** 2)]
    ])

    I_web = np.array([
        [(1 / 12) * m_w * (h ** 2 + t ** 2), 0, 0],
        [0, (1 / 12) * m_w * (L ** 2 + t ** 2), 0],
        [0, 0, (1 / 12) * m_w * (h ** 2 + L ** 2)]
    ])

    I_bot = np.array([
        [(1 / 12) * m_t * (s ** 2 + b ** 2), 0, 0],
        [0, (1 / 12) * m_t * (L ** 2 + b ** 2), 0],
        [0, 0, (1 / 12) * m_t * (s ** 2 + L ** 2)]
    ])

    displacement_top = np.array([0, 0, (d-s)/2])
    displacement_bot = np.array([0, 0, -(d-s)/2])
    translation_top = m_t * (
                np.dot(displacement_top, displacement_top) * np.identity(3) - np.outer(displacement_top, displacement_top))
    translation_bot = m_t * (
                np.dot(displacement_bot, displacement_bot) * np.identity(3) - np.outer(displacement_bot, displacement_bot))
    # https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor

    I = I_top + I_web + I_bot + translation_top + translation_bot

    return I

