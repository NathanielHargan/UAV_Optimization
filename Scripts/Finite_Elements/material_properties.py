# Material Properties (hopefully N, mm, MPa )
import math
Gravity = 9806.6499994  # mm/s^2

Materials = {
    # https://matweb.com/search/datasheet.aspx?matguid=4f19a42be94546b686bbf43f79c51b7d
    "Aluminum7075-T6": {
        "weight density": 0.0,  # [N/mm^3]
        "mass density": 2.81e-9,  # [N*s^2/mm^4] (N-s^2/mm)/mm^3 Mg/mm^3
        "elastic modulus": 71700,  # [MPa]
        "shear modulus": 26954.88721804511,  # [MPa]
        "poisson ratio": 0.33,  # [-]
        "shear yield strength": 503 * 0.577,  # [MPa]
        "tension yield strength": 503,  # [MPa]
        "ultimate strength": 572,  # [MPa]
        "fatigue strength coefficient": 1466,  # [MPa]
        "fatigue strength exponent": -0.143,  # [-]
        "fatigue strength cycles": 1e6  # [-]
    },
    "Carbon Fiber": {
        "weight density": 0.0,  # [N/mm^3]
        "mass density": 1.2179158e-9,  # [N*s^2/mm^4] (N-s^2/mm)/mm^3 Mg/mm^3
        "elastic modulus": 19000,  # [MPa]
        "shear modulus":  0.5 * 19000 / 1.3,  # [MPa]
        "poisson ratio": 0.3,  # calculate,  # [-]
        "shear yield strength": 345 * 0.577,  # [MPa]
        "tension yield strength": 345,  # [MPa]
        "ultimate strength": 345,  # [MPa]
        "fatigue strength coefficient": 1466,  # [MPa] # Change later
        "fatigue strength exponent": -0.143,  # [-] # Change later
        "fatigue strength cycles": 1e6,  # [-] # Change later
        "completely_reversed_stress_amplitude": 0
    }
}

# Calculate Fatigue Strength exponent and coefficient
Materials["Carbon Fiber"]["fatigue strength exponent"] = math.log(4/9) / math.log(Materials["Carbon Fiber"]["fatigue strength cycles"])
Materials["Carbon Fiber"]["fatigue strength coefficient"] = 0.9 * Materials["Carbon Fiber"]["ultimate strength"]

Materials["Carbon Fiber"]["completely_reversed_stress_amplitude"] = Materials["Carbon Fiber"]["fatigue strength coefficient"]  * (2 * Materials["Carbon Fiber"]["fatigue strength cycles"]) ** (Materials["Carbon Fiber"]["fatigue strength exponent"])

'''
# Tests 

print(Materials["Carbon Fiber"]["fatigue strength coefficient"])
print(Materials["Carbon Fiber"]["fatigue strength exponent"])
print(Materials["Carbon Fiber"]["ultimate strength"] * 0.9 * (2 * 0.5) ** Materials["Carbon Fiber"]["fatigue strength exponent"])
print(Materials["Carbon Fiber"]["ultimate strength"] * 0.9 * (2 * 5e8) ** Materials["Carbon Fiber"]["fatigue strength exponent"])
'''

Materials["Aluminum7075-T6"]["weight density"] = Materials["Aluminum7075-T6"]["mass density"] * Gravity
Materials["Carbon Fiber"]["weight density"] = Materials["Carbon Fiber"]["mass density"] * Gravity

def NU_Daniel_failure(normal_stress, shear_stress, Material):
    transverse_normal_tensile = Material["tension yield strength"]
    in_plane_shear_strength = Material["shear yield strength"]
    alpha = Material["elastic modulus"]/Material["shear modulus"]
    shear_failure_nu_d = (shear_stress/in_plane_shear_strength) ** 2 + normal_stress/in_plane_shear_strength * (2/alpha)
    tension_failure_nu_d = normal_stress/transverse_normal_tensile + (alpha/2) ** 2 * (shear_stress/transverse_normal_tensile) ** 2
    return shear_failure_nu_d, tension_failure_nu_d



#%%
