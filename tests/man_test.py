import sys
sys.path.append("../")

from thermoeos.mixture import Mixture
from thermoeos.equations import SRK
from thermoeos.equilibrium import find_equilibrium
# import matplotlib.pyplot as plt

# data = {
#     "name": "benzene-toluene",
#     "x": [0.089],
#     "y": [0.089],
#     "Tb": [353.25, 383.75],
#     "Tc": [562.1, 591.8],
#     "Pc": [48.9e5, 41.0e5],
#     "omega": [0.212, 0.263]
# }

data = {
    "name": "benzene-aniline",
    "x": [0.01],
    "y": [0.01],
    "Tb": [353.25, 457.25],
    "Tc": [562.1, 705.0],
    "Pc": [48.9e5, 56.3e5],
    "omega": [0.212, 0.382]
}

mix = Mixture(**data)

srk = SRK(mix)

# res = srk.solve(390, 101325)

res = find_equilibrium(srk, 400, 50000)

print(res[0])
print(f"T: {res[0].T}")
print(f"P: {res[0].P}")
print(res[1])
