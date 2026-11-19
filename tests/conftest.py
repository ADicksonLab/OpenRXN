"""Shared fixtures for the OpenRXN test suite.

All models are built here from the public API so individual tests
exercise one consistent, correct construction path.
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.connections import IsotropicConnection
from openrxn.model import Model


@pytest.fixture(autouse=True)
def _seed_rng():
    # The Gillespie propagator draws from the global np.random state;
    # seeding here makes every stochastic test reproducible.
    np.random.seed(0)


@pytest.fixture
def decay_model():
    """Single compartment, A --kf--> (nothing).

    Analytic solution: A(t) = A0 * exp(-kf * t).
    """
    A = Species("A")
    kf = 0.5 / unit.sec
    rxn = Reaction("decay", [A], [], [1], [], kf=kf)
    c = Compartment("main")
    c.add_rxn_to_compartment(rxn)
    flat = Model(compartments=[c]).flatten()
    return {"flat": flat, "species": {"A": A}, "kf": kf}


@pytest.fixture
def equilibrium_model():
    """Single compartment, A <=kf/kr=> B.

    Steady state: B/A = kf/kr;  A + B conserved.
    """
    A = Species("A")
    B = Species("B")
    kf = 1.0 / unit.sec
    kr = 0.25 / unit.sec
    rxn = Reaction("AB", [A], [B], [1], [1], kf=kf, kr=kr)
    c = Compartment("main")
    c.add_rxn_to_compartment(rxn)
    flat = Model(compartments=[c]).flatten()
    return {"flat": flat, "species": {"A": A, "B": B}, "kf": kf, "kr": kr}


@pytest.fixture
def two_box_model():
    """Two compartments connected by symmetric diffusion of species A.

    No reactions: A only moves between boxes, so total A is conserved and
    the two boxes equilibrate to equal amounts.
    """
    A = Species("A")
    k = 0.3 / unit.sec
    left = Compartment("left", volume=1.0 * unit.L)
    right = Compartment("right", volume=1.0 * unit.L)
    conn = IsotropicConnection({"A": k})
    left.connect(right, conn)
    right.connect(left, conn)
    flat = Model(compartments=[left, right]).flatten()
    return {"flat": flat, "species": {"A": A}, "k": k}
