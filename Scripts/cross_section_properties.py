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
        'second moment of area x': (math.pi / 4) * (r_outer ** 4 - r_inner ** 4),
        'second moment of area y': (math.pi / 4) * (r_outer ** 4 - r_inner ** 4),
        'second moment of area z': (math.pi / 2) * (r_outer ** 4 - r_inner ** 4),  # also polar moment of inertia
        'radius of gyration x': (1/2) * math.sqrt(r_outer ** 2 + r_inner ** 2),
        'radius of gyration y': (1/2) * math.sqrt(r_outer ** 2 + r_inner ** 2),
        'radius of gyration z': (math.sqrt(2)/2) * math.sqrt(r_outer ** 2 + r_inner ** 2),
        'plastic section modulus': (4/3) * (r_outer ** 3 - r_inner ** 3),
        'elastic section modulus': (math.pi/4) * (r_outer ** 4 - r_inner **4) / r_outer,
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
        'second moment of area x': (math.pi / 4) * r ** 4,
        'second moment of area y': (math.pi / 4) * r ** 4,
        'second moment of area z': (math.pi / 2) * r ** 4,  # also polar moment of inertia
        'radius of gyration x': r/2,
        'radius of gyration y': r/2,
        'radius of gyration z': (math.sqrt(2)/2) * r,
        'plastic section modulus': (4/3) * r ** 3,
        'elastic section modulus': (math.pi/4) * r ** 3,
        'torsional constant': (math.pi/2) * r ** 4,
        'transverse shear deflection constant y': 6/7,
        'transverse shear deflection constant z': 6/7
    }