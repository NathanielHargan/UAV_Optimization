# Material Properties (hopefully m-kg-s)
import pandas as pd

Materials = {
    # https://www.matweb.com/search/datasheet.aspx?matguid=b8d536e0b9b54bd7b69e4124d8f1d20a&ckck=1
    "Aluminum7075-T6": {
        "density": 2753.23377521,  # [kg/m^3]
        "elastic modulus": 6.89e7,  # [Pa]
        "shear modulus": 2.6e7,  # [Pa]
        "poisson ratio": 0.33,  # [-]
        "yield strength": 2.76e8,  # [Pa]
        "ultimate strength": 3.1e8,  # [Pa]
        "fatigue strength coefficient": 1.466e9,  # [Pa]
        "fatigue strength exponent": -0.134,  # [-]
        "fatigue strength cycles": 5e8  # [-]
    }
}



