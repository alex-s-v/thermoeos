import sys
sys.path.append("../")

from thermoeos.mixture import Mixture
from thermoeos.equations import SRK
from thermoeos.equilibrium import find_equilibrium

data = {
    "name": "benzene-toluene",
    "x": [0.5],
    "y": [0.5],
    "Tb": [353.25, 383.75],
    "Tc": [562.1, 591.8],
    "Pc": [48.9e5, 41.0e5],
    "omega": [0.212, 0.263]
}

mix = Mixture(**data)

srk = SRK(mix)

# res = srk.solve(390, 101325)

res = find_equilibrium(srk, 400, 50000)

print(res[0])
print(res[1])
