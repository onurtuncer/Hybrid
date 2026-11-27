# Hybrid

[![CI](https://github.com/onurtuncer/Hybrid/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Hybrid/actions/workflows/ci.yml?query=branch%3Amain)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Oxygen dissociation at fixed enthalpy and pressure

[tests/integration/test_oxygen_hp.py](tests/integration/test_oxygen_hp.py) runs NASA CEA on the [reference input deck](config/oxygen_hp.inp). The model considers neutral atomic and molecular oxygen in an adiabatic gas at constant pressure.

From the repository root in a Python 3.12 virtual environment:

```sh
python -m pip install -r requirements.txt
python -m pytest tests/integration -v
```

The test checks convergence, pressure, composition normalization, enthalpy conservation and the temperature bound for endothermic dissociation. Edit pressure (bar) or inlet temperature (K) in `config/oxygen_hp.inp` and rerun; add `-s` to print the equilibrium temperature and mole fractions.
