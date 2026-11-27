"""Integration test: oxygen dissociation at fixed enthalpy and pressure.

Runs NASA CEA on config/oxygen_hp.inp, a small CEA-style HP deck subset with
one O2 inlet, and checks the neutral O2/O equilibrium. Enthalpy follows the
normalization in NASA's Python sample:
https://github.com/nasa/cea/blob/main/source/bind/python/cea/samples/h2_o2.py
"""

import math
import re
from dataclasses import dataclass
from pathlib import Path

import cea
import numpy as np
import pytest

DECK = Path(__file__).resolve().parents[2] / "config" / "oxygen_hp.inp"
NUMBER = r"([0-9]+(?:\.[0-9]+)?)"
DECK_PATTERN = (
    r"problem\s+hp\s+p,bar=" + NUMBER
    + r"\s+reactants\s+name=O2\s+moles=" + NUMBER
    + r"\s+t,k=" + NUMBER + r"\s+end"
)


@dataclass(frozen=True)
class Deck:
    pressure: float  # bar
    inlet_moles: float
    inlet_temperature: float  # K


def parse_deck(text: str) -> Deck:
    fields = re.fullmatch(DECK_PATTERN, text.strip(), re.I)
    if fields is None:
        raise ValueError("Expected an HP deck with one O2 inlet.")
    deck = Deck(*map(float, fields.groups()))
    if min(deck.pressure, deck.inlet_moles, deck.inlet_temperature) <= 0:
        raise ValueError("Pressure, moles and temperature must be positive.")
    return deck


@pytest.fixture(scope="module")
def deck() -> Deck:
    return parse_deck(DECK.read_text())


@pytest.fixture(scope="module")
def equilibrium(deck):
    feed = cea.Mixture(["O2"])
    products = cea.Mixture(["O2", "O"])
    weights = feed.moles_to_weights(np.array([deck.inlet_moles]))
    enthalpy = feed.calc_property(cea.ENTHALPY, weights, deck.inlet_temperature)
    solver = cea.EqSolver(products, reactants=feed)
    state = cea.EqSolution(solver)
    solver.solve(state, cea.HP, enthalpy / cea.R, deck.pressure, weights)
    print(f"Equilibrium temperature: {state.T:.2f} K at {state.P:g} bar")
    for species, fraction in state.mole_fractions.items():
        print(f"{species:4s} {fraction:.8f}")
    return products, enthalpy, state


def test_converges(equilibrium):
    _, _, state = equilibrium
    assert state.converged


def test_holds_pressure(deck, equilibrium):
    _, _, state = equilibrium
    assert math.isclose(state.P, deck.pressure, rel_tol=1e-6)


def test_mole_fractions_are_normalized(equilibrium):
    _, _, state = equilibrium
    fractions = state.mole_fractions
    assert set(fractions) == {"O2", "O"}
    assert all(x >= 0 for x in fractions.values())
    assert math.isclose(sum(fractions.values()), 1, abs_tol=1e-6)


def test_endothermic_dissociation_cools_gas(deck, equilibrium):
    _, _, state = equilibrium
    assert 0 < state.T <= deck.inlet_temperature + 1e-4


def test_conserves_enthalpy(equilibrium):
    products, inlet_enthalpy, state = equilibrium
    fractions = state.mole_fractions
    moles = np.array([fractions[name] for name in products.species_names])
    weights = products.moles_to_weights(moles)
    outlet_enthalpy = products.calc_property(cea.ENTHALPY, weights, state.T)
    assert math.isclose(outlet_enthalpy, inlet_enthalpy, rel_tol=1e-6)


@pytest.mark.parametrize(
    "text",
    [
        "problem tp p,bar=1 reactants name=O2 moles=1 t,k=3500 end",
        "problem hp p,bar=0 reactants name=O2 moles=1 t,k=3500 end",
    ],
)
def test_rejects_unsupported_decks(text):
    with pytest.raises(ValueError):
        parse_deck(text)
