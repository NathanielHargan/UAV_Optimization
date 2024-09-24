# Material Properties (hopefully N, mm, MPa )
Gravity = 9806.6499994  # mm/s^2

Materials = {
    # https://matweb.com/search/datasheet.aspx?matguid=4f19a42be94546b686bbf43f79c51b7d
    "Aluminum7075-T6": {
        "weight density": 26.6e-6,  # [N/mm^3]
        "mass density": 0.0,  # [N*s^2/mm^4] #(N-s^2/mm)/mm^3
        "elastic modulus": 71700,  # [MPa]
        "shear modulus": 26900,  # [MPa]
        "poisson ratio": 0.33,  # [-]
        "yield strength": 503,  # [MPa]
        "ultimate strength": 572,  # [MPa]
        "fatigue strength coefficient": 1466,  # [MPa]
        "fatigue strength exponent": -0.143,  # [-]
        "fatigue strength cycles": 5e8  # [-]
    }
}

Materials["Aluminum7075-T6"]["mass density"] = Materials["Aluminum7075-T6"]["weight density"] / Gravity


