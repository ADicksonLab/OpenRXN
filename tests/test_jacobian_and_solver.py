"""C2 gates:
1. The analytic Jacobian (_jac) must match a finite-difference Jacobian of the
   RHS on random states (correctness of the Jacobian).
2. The default stiff solver (Radau + analytic Jacobian) must agree with the
   explicit RK45 on the analytic models (the solver swap doesn't change answers,
   only cost/accuracy).
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem


def _fd_jacobian(sys, Q, eps_rel=1e-6):
    """Central finite-difference Jacobian of the RHS at Q."""
    n = Q.size
    J = np.zeros((n, n))
    f0 = sys._dQ_dt(0.0, Q)  # noqa: F841 (kept for clarity / sanity)
    for k in range(n):
        h = eps_rel * max(abs(Q[k]), 1.0)
        Qp = Q.copy(); Qp[k] += h
        Qm = Q.copy(); Qm[k] -= h
        J[:, k] = (sys._dQ_dt(0.0, Qp) - sys._dQ_dt(0.0, Qm)) / (2 * h)
    return J


def _assert_jac_matches_fd(sys, n_trials=10, scale=1e2, seed=0):
    rng = np.random.default_rng(seed)
    for _ in range(n_trials):
        Q = rng.random(sys.state.size) * scale
        J_analytic = sys._jac(0.0, Q)
        J_fd = _fd_jacobian(sys, Q)
        assert np.allclose(J_analytic, J_fd, rtol=1e-5, atol=1e-6)


def test_jacobian_decay(decay_model):
    """Order-1 Jacobian matches finite differences."""
    _assert_jac_matches_fd(ODESystem(decay_model["flat"]))


def test_jacobian_equilibrium(equilibrium_model):
    """Reversible A<=>B Jacobian matches finite differences."""
    _assert_jac_matches_fd(ODESystem(equilibrium_model["flat"]))


def test_jacobian_two_box(two_box_model):
    """Diffusion (linear) Jacobian matches finite differences."""
    _assert_jac_matches_fd(ODESystem(two_box_model["flat"]))


def test_jacobian_bimolecular():
    """Order-2 A+B->C (and reverse) Jacobian matches finite differences."""
    A, B, C = Species("A"), Species("B"), Species("C")
    r = Reaction("ABC", [A, B], [C], [1, 1], [1], kf=0.01 / unit.sec, kr=0.2 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(r)
    _assert_jac_matches_fd(ODESystem(Model(compartments=[c]).flatten()))


def test_jacobian_dimerization():
    """2A->C (repeated reactant) Jacobian matches finite differences."""
    A, C = Species("A"), Species("C")
    r = Reaction("AAC", [A], [C], [2], [1], kf=0.01 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(r)
    _assert_jac_matches_fd(ODESystem(Model(compartments=[c]).flatten()))


def test_radau_matches_rk45_decay(decay_model):
    """Default Radau and explicit RK45 give the same decay trajectory."""
    sysd = ODESystem(decay_model["flat"])
    i = sysd.state.index["main"]["A"]
    sysd.set_q([i], 100.0)
    radau = sysd.run(8.0)  # default method='Radau'

    sysr = ODESystem(decay_model["flat"])
    sysr.set_q([i], 100.0)
    rk45 = sysr.run(8.0, method="RK45")

    k = decay_model["kf"].to(1 / unit.sec).magnitude
    truth = 100.0 * np.exp(-k * 8.0)
    assert radau.y[i][-1] == pytest.approx(truth, rel=1e-5)
    assert rk45.y[i][-1] == pytest.approx(truth, rel=1e-5)
    assert radau.y[i][-1] == pytest.approx(rk45.y[i][-1], rel=1e-4)


def test_radau_matches_rk45_equilibrium(equilibrium_model):
    """Default Radau and explicit RK45 reach the same A<=>B steady state."""
    iA_name, iB_name = "A", "B"
    sysd = ODESystem(equilibrium_model["flat"])
    iA = sysd.state.index["main"][iA_name]
    iB = sysd.state.index["main"][iB_name]
    sysd.set_q([iA], 100.0); sysd.set_q([iB], 0.0)
    radau = sysd.run(50.0)

    sysr = ODESystem(equilibrium_model["flat"])
    sysr.set_q([iA], 100.0); sysr.set_q([iB], 0.0)
    rk45 = sysr.run(50.0, method="RK45")

    assert radau.y[iB][-1] / radau.y[iA][-1] == pytest.approx(
        rk45.y[iB][-1] / rk45.y[iA][-1], rel=1e-4
    )
