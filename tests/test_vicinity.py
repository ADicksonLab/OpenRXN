"""Vicinity (ligand rebinding) rate law in ODESystem:
1. RHS matches the closed-form own- and shared-vicinity equations
   (Vauquelin 2016, Table 1 eq. 7; own: (kon D R - koff DR) / (1 + kon R / k_drift)).
2. A large k_drift recovers plain mass action.
3. The analytic Jacobian matches finite differences.
4. GillespieSystem raises NotImplementedError.
5. Invalid Vicinity / tagged Reaction inputs are rejected.
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction, Vicinity
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem
from openrxn.systems.GillespieSystem import GillespieSystem

NA = 6.022e23
KON = np.array([1e8, 1e6])
KOFF = 1e-9 * KON
KE = 1e-3


def _build(k_drift=1e-2, shared=False, release=False, system=ODESystem):
    """D + R_i <-> DR_i (i = 1, 2) in 1 L, plus export of D; release adds DR_i -> D into the vicinity."""
    D, R1, R2, DR1, DR2 = (Species(n) for n in ['D', 'R1', 'R2', 'DR1', 'DR2'])
    vic = [Vicinity('v', D, k_drift / unit.sec)] * 2 if shared else \
          [Vicinity(f'v{i}', D, k_drift / unit.sec) for i in (1, 2)]
    rxns = [Reaction('export', [D], [], [1], [], kf=1e-2 / unit.sec)]
    for i, (R, DR) in enumerate([(R1, DR1), (R2, DR2)]):
        rxns.append(Reaction(f'bind{i}', [D, R], [DR], [1, 1], [1], kf=KON[i] / (unit.molar * unit.sec),
                             kr=KOFF[i] / unit.sec, vicinity=vic[i]))
        if release:
            rxns.append(Reaction(f'rel{i}', [DR], [D], [1], [1], kf=KE / unit.sec, vicinity=vic[i]))
    c = Compartment('main', volume=1 * unit.L)
    c.add_rxns_to_compartment(rxns)
    sys = system(Model(compartments=[c]).flatten())
    return sys, [list(sys.state.species).index(n) for n in ['D', 'R1', 'R2', 'DR1', 'DR2']]


def _random_states(n_trials=20, seed=0):
    rng = np.random.default_rng(seed)
    for _ in range(n_trials):
        yield rng.uniform(1e-10, 2e-7), rng.uniform(0, 2e-9, 2), rng.uniform(0, 2e-9, 2)   # D, R, DR (M)


def _rhs_molar(sys, idx, D, R, DR):
    Q = np.zeros(sys.state.size)
    Q[idx] = np.concatenate([[D], R, DR]) * NA
    return sys._dQ_dt(0.0, Q)[idx] / NA


@pytest.mark.parametrize("release", [False, True])
def test_own_vicinity_closed_form(release):
    sys, idx = _build(release=release)
    ke = KE if release else 0.0
    for D, R, DR in _random_states():
        v = (KON * D * R - (KOFF + ke) * DR) / (1 + KON * R / 1e-2)
        dR = -v - ke * DR
        expected = np.concatenate([[-1e-2 * D - v.sum()], dR, v])
        assert np.allclose(_rhs_molar(sys, idx, D, R, DR), expected, rtol=1e-9, atol=1e-30)


def test_shared_vicinity_closed_form():
    sys, idx = _build(shared=True)
    kd = 1e-2
    for D, R, DR in _random_states():
        den = kd + (KON * R).sum()
        hop = (KON[0] * R[0] * KOFF[1] * DR[1] - KON[1] * R[1] * KOFF[0] * DR[0]) / den
        v = kd * (KON * D * R - KOFF * DR) / den + np.array([hop, -hop])
        expected = np.concatenate([[-1e-2 * D - v.sum()], -v, v])
        assert np.allclose(_rhs_molar(sys, idx, D, R, DR), expected, rtol=1e-9, atol=1e-30)


@pytest.mark.parametrize("shared", [False, True])
def test_large_k_drift_is_mass_action(shared):
    sys, idx = _build(k_drift=1e12, shared=shared)
    for D, R, DR in _random_states():
        v = KON * D * R - KOFF * DR
        expected = np.concatenate([[-1e-2 * D - v.sum()], -v, v])
        assert np.allclose(_rhs_molar(sys, idx, D, R, DR), expected, rtol=1e-6, atol=1e-30)


@pytest.mark.parametrize("shared", [False, True])
def test_jacobian_matches_fd(shared):
    sys, idx = _build(shared=shared, release=True)
    for D, R, DR in _random_states(n_trials=5):
        Q = np.zeros(sys.state.size)
        Q[idx] = np.concatenate([[D], R, DR]) * NA
        J_fd = np.zeros((Q.size, Q.size))
        for k in range(Q.size):
            h = 1e-6 * max(abs(Q[k]), 1.0)
            Qp, Qm = Q.copy(), Q.copy()
            Qp[k] += h
            Qm[k] -= h
            J_fd[:, k] = (sys._dQ_dt(0.0, Qp) - sys._dQ_dt(0.0, Qm)) / (2 * h)
        assert np.allclose(sys._jac(0.0, Q), J_fd, rtol=1e-5, atol=1e-12)


def test_gillespie_not_implemented():
    with pytest.raises(NotImplementedError):
        _build(system=GillespieSystem)


def test_invalid_inputs():
    D, R, DR = Species('D'), Species('R'), Species('DR')
    kon = 1e8 / (unit.molar * unit.sec)
    with pytest.raises(TypeError):
        Vicinity('v', 'D', 1e-2 / unit.sec)                       # ligand is not a Species
    with pytest.raises(ValueError):
        Vicinity('v', D, 0 / unit.sec)                            # k_drift must be positive
    v = Vicinity('v', D, 1e-2 / unit.sec)
    assert repr(v) == 'Vicinity(v)'
    with pytest.raises(ValueError):
        Reaction('deg', [R], [], [1], [], kf=1 / unit.sec, vicinity=v)                       # ligand absent
    with pytest.raises(ValueError):
        Reaction('b2', [D, R], [DR], [2, 1], [1], kf=kon / unit.molar, vicinity=v)          # stoichiometry 2
    with pytest.raises(ValueError):
        Reaction('both', [D, R], [D, DR], [1, 1], [1, 1], kf=kon, vicinity=v)               # on both sides
