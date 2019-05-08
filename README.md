# About

ThermoEoS is the python package for engineers and scientis which can be used 
for calculation of equlibrium properties of the multicomponent mixtures.

# Usage examples

```python
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

mix = Mixture(**data) # Creating the mixture
srk = SRK(mix)        # Optimizing EoS for this mixture

# Solving EoS
new_mix, issuccessful = srk.solve(T=400, P=1013250)

print(mew_mix)
print(issuccessful)

# Calculating equilibrium pressure and molar fraction of the components
# in the vapor phase for T = 400 K 
eqil_mix, issuccessful = find_equilibrium(srk, T=400, P=50000)

print(eqil_mix)
print(issuccessful)
```

# Notes

The package currently in an early development stage and major changes can occur.

You can use this package along with [Caleb Bell's thermo package](https://github.com/CalebBell/thermo/),
for better user experience.
