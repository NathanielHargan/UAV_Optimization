import math
import numpy as np

# returns diction annulus beam cross section properties
def cross_section_annulus(cross_section_parameters):
    r_outer = cross_section_parameters[0]/2
    r_inner = cross_section_parameters[1]/2
    m = r_inner/r_outer
    # https://structx.com/Shape_Formulas_014.html
    return {
        'area': math.pi * (r_outer ** 2 - r_inner ** 2),
        'perimeter': math.pi * 2 * (r_outer + r_inner),
        'outer perimeter': math.pi * 2 * r_outer,
        'inner perimeter': math.pi * 2 * r_inner,
        'second moment of area x': (math.pi / 2) * (r_outer ** 4 - r_inner ** 4),  # also polar moment of inertia
        'second moment of area y': (math.pi / 4) * (r_outer ** 4 - r_inner ** 4),
        'second moment of area z': (math.pi / 4) * (r_outer ** 4 - r_inner ** 4),
        'radius of gyration x': (math.sqrt(2)/2) * math.sqrt(r_outer ** 2 + r_inner ** 2),
        'radius of gyration y': (1/2) * math.sqrt(r_outer ** 2 + r_inner ** 2),
        'radius of gyration z': (1/2) * math.sqrt(r_outer ** 2 + r_inner ** 2),
        'plastic section modulus': (4/3) * (r_outer ** 3 - r_inner ** 3),
        'elastic section modulus': (math.pi/4) * (r_outer ** 4 - r_inner ** 4) / r_outer,
        'torsional constant': (math.pi/2) * (r_outer ** 4 - r_inner ** 4),
        # I calculated this for poisson = 0
        'transverse shear deflection constant y': (6 * (1 + m ** 2) ** 2) / ((7 * (1 + m ** 2) ** 2) + (20 * m ** 2)),
        'transverse shear deflection constant z': (6 * (1 + m ** 2) ** 2) / ((7 * (1 + m ** 2) ** 2) + (20 * m ** 2)),
        'top fiber y': r_outer,
        'top fiber z': r_outer,
        'stress points': [],
        'circumscribed': r_outer
    }


def cross_section_circle(cross_section_parameters):
    r = cross_section_parameters[0]/2

    return {
        'area': math.pi * (r ** 2),
        'perimeter': math.pi * 2 * r,
        'second moment of area x': (math.pi / 2) * r ** 4,  # also polar moment of inertia
        'second moment of area y': (math.pi / 4) * r ** 4,
        'second moment of area z': (math.pi / 4) * r ** 4,
        'radius of gyration x': (math.sqrt(2)/2) * r,
        'radius of gyration y': r/2,
        'radius of gyration z': r/2,
        'plastic section modulus': (4/3) * r ** 3,
        'elastic section modulus': (math.pi/4) * r ** 3,
        'torsional constant': (math.pi/2) * r ** 4,
        'transverse shear deflection constant y': 6/7,
        'transverse shear deflection constant z': 6/7,
        'stress points': [],
        'circumscribed': r
    }


def cross_section_rectangle(cross_section_parameters):
    # https://structx.com/Shape_Formulas_024.html
    h = cross_section_parameters[0]  # height
    b = cross_section_parameters[1]  # base
    # h in y dir b in z dir
    return {
        'area': h*b,
        'perimeter': 2 * (h + b),
        'second moment of area x': h * b * ((h ** 2) + (b ** 2)) / 12,  # also polar moment of inertia
        'second moment of area y': ((b ** 3) * h) / 12,  # I1
        'second moment of area z': (b * (h ** 3)) / 12,  # I2
        'radius of gyration x': math.sqrt((h ** 2 + b ** 2) / (2 * math.sqrt(3))),
        'radius of gyration y': b / (2 * math.sqrt(3)),
        'radius of gyration z': h / (2 * math.sqrt(3)),
        'plastic section modulus x': h ** 2 * b / 4,
        'plastic section modulus y': h * b ** 2 / 4,
        'elastic section modulus x': (h ** 2 * b) / 6,
        'elastic section modulus y': (h * b ** 2) / 6,
        'torsional constant': ((h ** 3) * b) * ((1/3) - (0.21*h/b) * (1 - (h**4 / (12 * b ** 4)))),
        'transverse shear deflection constant y': 5/6,
        'transverse shear deflection constant z': 5/6,
        'stress points': np.array([[0,h/2],[b/2,h/2],[b/2,0],[b/2,-h/2],[0,-h/2],[-b/2,-h/2],[-b/2,0],[-b/2,h/2]]),
        'circumscribed': np.linalg.norm([h/2,b/2])
        # 'stress points': np.array([[h/2,0],[0,b/2],[-h/2,0],[0,-b/2]])
    }

def cross_section_hexagon(cross_section_parameters):
    # https://structx.com/Shape_Formulas_036.html
    r = cross_section_parameters[0]/2  # radius
    f = r * math.sqrt(3)  # face to face distance
    a = (3 * math.sqrt(3) * r ** 2) / 2  # area
    iz = 0.0601 * f ** 4
    return {
        'area': a,
        'perimeter': 6 * r,
        'second moment of area x': 2 * iz,  # also polar moment of inertia
        'second moment of area y': iz,  # I1
        'second moment of area z': iz,  # I2
        'radius of gyration x': math.sqrt(2 * iz/a),
        'radius of gyration y': math.sqrt(iz/a),
        'radius of gyration z': math.sqrt(iz/a),
        'elastic section modulus': iz/r,
        'torsional constant': 0.1154 * f ** 4,
        'transverse shear deflection constant y': 0.00001,
        'transverse shear deflection constant z': 0.00001,
        'stress points': np.array([[r,0], [f/2,r/2], [f/2,-r/2], [0,-r], [-f/2,r/2], [-f/2,-r/2]]),
        'circumscribed': r
    }


# https://structx.com/Shape_Formulas_013.html x=z
def cross_section_i_beam(cross_section_parameters):
    b = cross_section_parameters[0]
    d = cross_section_parameters[1]
    s = cross_section_parameters[2]
    t = cross_section_parameters[3]

    h = d - 2*s
    cz = b/2
    cy = d/2
    # Area
    area = b * d - h * (b - t)

    # Second moment strong axis (z-axis, horizontal bending)
    Iz = (b * d ** 3) / 12.0 - ((b - t) * h ** 3) / 12.0

    # Second moment weak axis (y-axis, vertical bending)
    Iy_flange = 2 * (s * b ** 3) / 12.0
    Iy_web = (h * t ** 3) / 12.0
    Iy = Iy_flange + Iy_web

    # Polar second moment (centroidal)
    J_polar = Iy ** 2 + Iz

    # Saint-Venant torsional constant
    J_torsion = (1.0 / 3.0) * (2 * b * s ** 3 + h * t ** 3)

    # Radii of gyration
    kx = math.sqrt(J_polar / area)
    ky = math.sqrt(Iy / area)
    kz = math.sqrt(Iz / area)

    # Elastic section modulus
    Sy = Iy / (b / 2.0)
    Sz = Iz / (d / 2.0)

    # Plastic section moduli
    Zp_z = 2 * (b * s * (d / 2.0 - s / 2.0) + t * (h / 2.0) * (h / 4.0))
    Zp_y = (s * b ** 2) / 2.0 + (h * t ** 2) / 4.0

    # Shear deflection constants
    m = (2*b*s)/(h*t)
    n = b / h
    kappa_z = 10 * (1 + 3*m) ** 2 / ((12+72*m+150*m**2+90*m**3) + (30*n**2)*(m+m**2))
    kappa_y = 10 * (1 + 3*m) ** 2 / ((12+72*m+150*m**2+90*m**3) + (30*n**2)*(m+m**2))

    return {
        'area': area,
        'second moment of area y': Iy,
        'second moment of area z': Iz,
        'second moment of area x': J_polar,
        'polar moment of inertia': J_polar,
        'torsional constant': J_torsion,
        'radius of gyration x': kx,
        'radius of gyration y': ky,
        'radius of gyration z': kz,
        'elastic section modulus y': Sy,
        'elastic section modulus z': Sz,
        'plastic section modulus y': Zp_y,
        'plastic section modulus z': Zp_z,
        'transverse shear deflection constant y': kappa_y,
        'transverse shear deflection constant z': kappa_z,
        'top fiber y': cy,
        'top fiber z': cz,
        'circumscribed': np.sqrt(cy ** 2 + cz ** 2)
    }
