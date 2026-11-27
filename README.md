# Hybrid

[01 · Oxygen dissociation at fixed enthalpy and pressure](notebooks/01_oxygen_hp.ipynb) runs NASA CEA on the [reference input deck](config/oxygen_hp.inp). The model considers neutral atomic and molecular oxygen in an adiabatic gas at constant pressure.

## Run the study

From the repository root in a Python 3.12 virtual environment:

```sh
python -m pip install -r requirements.txt
jupyter nbconvert --execute --to html notebooks/01_oxygen_hp.ipynb --output-dir=build
```

Edit pressure (bar) or inlet temperature (K) in `config/oxygen_hp.inp` and rerun. The notebook reports equilibrium temperature and mole fractions and checks convergence, pressure, composition normalization and the temperature bound for endothermic dissociation.
