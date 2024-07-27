import math


# returns diction annulus beam cross section properties
def cross_section_annulus(cross_section_parameters):
    r_outer = cross_section_parameters[0]
    r_inner = cross_section_parameters[1]
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
        'transverse shear deflection constant y': (6 * (1 + m) ** 2) / ((7 * (1 + m) ** 2) + (20 * m ** 2)),
        'transverse shear deflection constant z': (6 * (1 + m) ** 2) / ((7 * (1 + m) ** 2) + (20 * m ** 2))
    }


def cross_section_circle(cross_section_parameters):
    r = cross_section_parameters[0]

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
        'transverse shear deflection constant z': 6/7
    }


def cross_section_rectangle(cross_section_parameters):
    # https://structx.com/Shape_Formulas_024.html
    h = cross_section_parameters[0]  # height
    b = cross_section_parameters[1]  # base

    return {
        'area': h*b,
        'perimeter': 2 * (h + b),
        'second moment of area x': h * b * ((h ** 2) + (b ** 2)) / 12,  # also polar moment of inertia
        'second moment of area y': ((b ** 3) * h) / 12,  # I1
        'second moment of area z': (b * (h ** 3)) / 12,  # I2
        'radius of gyration x': h / (2 * math.sqrt(3)),
        'radius of gyration y': b / (2 * math.sqrt(3)),
        'radius of gyration z': math.sqrt((h ** 2 + b ** 2) / (2 * math.sqrt(3))),
        'plastic section modulus x': h ** 2 * b / 4,
        'plastic section modulus y': h * b ** 2 / 4,
        'elastic section modulus x': (h ** 2 * b) / 6,
        'elastic section modulus y': (h * b ** 2) / 6,
        'torsional constant': ((h ** 3) * b) * ((1/3) - (0.21*h/b) * (1 - (h**4 / (12 * b ** 4)))),
        'transverse shear deflection constant y': 5/6,
        'transverse shear deflection constant z': 5/6
    }

def cross_section_hexagon(cross_section_parameters):
    # https://structx.com/Shape_Formulas_036.html
    r = cross_section_parameters[0]  # radius
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
        'transverse shear deflection constant y': 0,
        'transverse shear deflection constant z': 0
    }


def cross_section_octagon(cross_section_parameters):
    # https://structx.com/Shape_Formulas_037.html
    r = cross_section_parameters[0]  # radius
    return {
        'area': 4 * r ** 2 * math.sqrt(2) / 2,
        'perimeter': 8 * r * math.sqrt(2 - math.sqrt(2)),
        'second moment of area x': 1.2762 * r ** 4,  # also polar moment of inertia
        'second moment of area y': 0.6381 * r ** 4,  # I1
        'second moment of area z': 0.6381 * r ** 4,  # I2
        'radius of gyration x': 0.672 * r,
        'radius of gyration y': 0.475 * r,
        'radius of gyration z': 0.475 * r,
        'elastic section modulus': 0.6381 * r ** 3,
        'transverse shear deflection constant y': 0,
        'transverse shear deflection constant z': 0
    }