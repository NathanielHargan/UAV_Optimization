# Material Properties (hopefully N, mm, MPa )
Gravity = 9806.6499994  # mm/s^2

Materials = {
    # https://matweb.com/search/datasheet.aspx?matguid=4f19a42be94546b686bbf43f79c51b7d
    "Aluminum7075-T6": {
        "weight density": 2.81e-5,  # [N/mm^3]
        "mass density": 0,  # [N-s^2/mm^4]
        "elastic modulus": 71.7,  # [MPa]
        "shear modulus": 26.9,  # [MPa]
        "poisson ratio": 0.33,  # [-]
        "yield strength": 503,  # [MPa]
        "ultimate strength": 572,  # [MPa]
        "fatigue strength coefficient": 1466,  # [MPa]
        "fatigue strength exponent": -0.143,  # [-]
        "fatigue strength cycles": 5e8  # [-]
    }
}

Materials["Aluminum7075-T6"]["mass density"] = Materials["Aluminum7075-T6"]["weight density"] / Gravity



