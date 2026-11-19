"""C1 parity gate: the vectorized RHS (ODESystem._dQ_dt) must reproduce the
original per-state-variable RHS (ODESystem._dQ_dt_original) on arbitrary state
vectors. The two differ only in floating-point summation order, so we compare
with a tight tolerance.

Covers reaction orders 0 (birth), 1 (decay/equilibrium), 2 (bimolecular),
and the inter-compartment diffusion path (two_box).
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem


def _assert_rhs_parity(sys, n_trials=25, scale=1e3, seed=0):
    rng = np.random.default_rng(seed)
    for _ in range(n_trials):
        Q = rng.random(sys.state.size) * scale
        new = sys._dQ_dt(0.0, Q)
        original = sys._dQ_dt_original(0.0, Q)
        assert np.allclose(new, original, rtol=1e-10, atol=1e-12)


def test_rhs_parity_decay(decay_model):
    """Order-1 decay: vectorized RHS matches legacy."""
    _assert_rhs_parity(ODESystem(decay_model["flat"]))


def test_rhs_parity_equilibrium(equilibrium_model):
    """Order-1 reversible A<=>B: vectorized RHS matches legacy."""
    _assert_rhs_parity(ODESystem(equilibrium_model["flat"]))


def test_rhs_parity_two_box(two_box_model):
    """Diffusion (connection) terms: vectorized RHS matches legacy."""
    _assert_rhs_parity(ODESystem(two_box_model["flat"]))


def test_rhs_parity_bimolecular():
    """Order-2 A+B->C (reverse C->A+B) exercises the product-of-two path."""
    A, B, C = Species("A"), Species("B"), Species("C")
    r = Reaction("ABC", [A, B], [C], [1, 1], [1], kf=0.01 / unit.sec, kr=0.2 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(r)
    _assert_rhs_parity(ODESystem(Model(compartments=[c]).flatten()))


def test_rhs_parity_birth_death():
    """Order-0 birth (reverse of A->nothing) exercises the constant-term path."""
    A = Species("A")
    r = Reaction("bd", [A], [], [1], [], kf=0.1 / unit.sec, kr=1.0 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(r)
    _assert_rhs_parity(ODESystem(Model(compartments=[c]).flatten()))


def test_rhs_parity_dimerization():
    """Order-2 with repeated reactant 2A->C exercises prod(Q[[i,i]])."""
    A, C = Species("A"), Species("C")
    r = Reaction("AAC", [A], [C], [2], [1], kf=0.01 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(r)
    _assert_rhs_parity(ODESystem(Model(compartments=[c]).flatten()))
