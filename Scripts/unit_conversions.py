from Scripts.Finite_Elements.material_properties import Gravity
import math

# Length
mm_to_m = 1000

# Mass
kg_to_Mg = 10 ** -3
g_to_Mg = 10 ** -6
Mg_to_lb = 2.20462*1000


# Power
mW_to_kW = 10 ** -6
kW_to_mW = 10 ** 6

# Time
hrs_to_s = 3600
s_to_hrs = 1/3600

# Thrust
N_to_kgf = 1000/(Gravity)

rad_per_s_to_Hz = 1/(2*math.pi)
Hz_to_rad_per_s = 2*math.pi

rpm_to_rad_per_s = math.pi / 30

