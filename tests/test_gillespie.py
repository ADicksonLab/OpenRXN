"""Stochastic (Gillespie) validation: immigration-death mean vs analytic
steady state, exact count-conservation invariants, set_q unit handling,
volume scaling of rates, and transport through connections. RNG is seeded by
the autouse fixture in conftest, so results are reproducible.
Also checks bimolecular rates against the ODE (Avogadro factor), DivByV
transport, runs that reach a state where nothing can fire, and that the first
event is not fired at t0.
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.connections import DivByVConnection
from openrxn.model import Model
from openrxn.systems.GillespieSystem import GillespieSystem
from openrxn.systems.ODESystem import ODESystem
from openrxn.systems.reporters import AllReporter


def test_birth_death_mean_matches_ode():
    """Immigration-death A reaches stochastic mean kr/kf = 10 (averaged over runs)."""
    A = Species("A")
    rxn = Reaction("bd", [A], [], [1], [], kf=0.1 / unit.sec, kr=1.0 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(rxn)
    flat = Model(compartments=[c]).flatten()

    finals = []
    for _ in range(100):
        g = GillespieSystem(flat)
        i = g.state.index["main"]["A"]
        g.set_q([i], 0)
        g.run(60.0)
        finals.append(g.state.q_val[i])

    assert np.mean(finals) == pytest.approx(10.0, abs=2.0)


def test_gillespie_conserves_counts():
    """A<=>B is closed, so A + B is conserved exactly through a stochastic run."""
    A, B = Species("A"), Species("B")
    rxn = Reaction("AB", [A], [B], [1], [1], kf=1.0 / unit.sec, kr=0.5 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(rxn)
    flat = Model(compartments=[c]).flatten()

    g = GillespieSystem(flat)
    iA = g.state.index["main"]["A"]
    iB = g.state.index["main"]["B"]
    g.set_q([iA], 100)
    g.set_q([iB], 0)
    g.run(5.0)
    assert g.state.q_val[iA] + g.state.q_val[iB] == pytest.approx(100.0)


def test_gillespie_dimerization_parity():
    """2A --> C removes two A and adds one C per event, so A + 2C = A0 holds."""
    A, C = Species("A"), Species("C")
    rxn = Reaction("AAC", [A], [C], [2], [1], kf=0.01 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(rxn)
    flat = Model(compartments=[c]).flatten()

    g = GillespieSystem(flat)
    iA = g.state.index["main"]["A"]
    iC = g.state.index["main"]["C"]
    g.set_q([iA], 100)
    g.set_q([iC], 0)
    g.run(10.0)
    assert g.state.q_val[iA] + 2 * g.state.q_val[iC] == pytest.approx(100.0)


def _flat(rxns, volume=None):
    """Flat model with one compartment 'main' holding rxns."""
    c = Compartment("main", volume=volume)
    c.add_rxns_to_compartment(rxns)
    return Model(compartments=[c]).flatten()


def test_set_q_accepts_dimensionless_quantity(decay_model):
    """set_q strips a dimensionless Pint quantity down to its count."""
    g = GillespieSystem(decay_model["flat"])
    i = g.state.index["main"]["A"]
    g.set_q([i], 7 * unit.dimensionless)
    assert g.state.q_val[i] == 7.0


def test_set_q_rejects_dimensional_quantity(decay_model):
    """Gillespie counts cannot be set from mol or mol/L; set_q raises ValueError."""
    g = GillespieSystem(decay_model["flat"])
    i = g.state.index["main"]["A"]
    for q in (1.0 * unit.mol, 1.0 * unit.molar):
        with pytest.raises(ValueError, match="dimensionless"):
            g.set_q([i], q)
    assert g.state.q_val[i] == 0.0


def test_unimolecular_rate_ignores_volume():
    """A first-order rate is not volume-scaled: same seed, same trajectory with or without a volume."""
    A, B = Species("A"), Species("B")
    rxn = Reaction("AB", [A], [B], [1], [1], kf=0.5 / unit.sec, kr=0.2 / unit.sec)
    finals = []
    for vol in (None, 1e-15 * unit.L):
        np.random.seed(1)
        g = GillespieSystem(_flat([rxn], volume=vol))
        g.set_q([g.state.index["main"]["A"]], 50)
        g.run(2.0)
        finals.append(g.state.q_val[g.state.index["main"]["A"]])
    assert finals[0] == finals[1]
    assert 0 < finals[0] < 50


def test_reports_land_on_checkpoints(decay_model):
    """Each propagate stops exactly at its checkpoint, so reports are at t = 0, 1, 2, 3 and counts only fall."""
    g = GillespieSystem(decay_model["flat"])
    i = g.state.index["main"]["A"]
    g.set_q([i], 50)
    rep = AllReporter(freq=1)
    g.add_reporter(rep)
    out = g.run(3.0)
    assert [r["t"] for r in rep.reports()] == [0.0, 1.0, 2.0, 3.0]
    assert out["final_t"] == 3.0
    counts = [r["report"][i] for r in rep.reports()]
    assert counts[0] == 50
    assert all(a >= b for a, b in zip(counts, counts[1:]))


def test_isotropic_transport_relaxation(two_box_model):
    """All A starts left; mean left count follows N/2 (1 + exp(-2kt)) and each run conserves N."""
    N, t, k = 100, 1.0, two_box_model["k"].magnitude
    lefts = []
    for _ in range(30):
        g = GillespieSystem(two_box_model["flat"])
        iL = g.state.index["left"]["A"]
        iR = g.state.index["right"]["A"]
        g.set_q([iL], N)
        g.run(t)
        assert g.state.q_val[iL] + g.state.q_val[iR] == N
        lefts.append(g.state.q_val[iL])
    assert np.mean(lefts) == pytest.approx(N / 2 * (1 + np.exp(-2 * k * t)), abs=3.0)


def test_divbyv_transport_relaxation():
    """A DivByV connection hops at k = kV / V per molecule, like the ODE system."""
    left = Compartment("left", volume=1000 * unit.nm**3)
    right = Compartment("right", volume=1000 * unit.nm**3)
    conn = DivByVConnection({"A": 100 * unit.nm**3 / unit.sec})
    left.connect(right, conn)
    right.connect(left, conn)
    flat = Model(compartments=[left, right]).flatten()
    N, t, k = 100, 2.5, 0.1
    lefts = []
    for _ in range(30):
        g = GillespieSystem(flat)
        g.set_q([g.state.index["left"]["A"]], N)
        g.run(t)
        lefts.append(g.state.q_val[g.state.index["left"]["A"]])
    assert np.mean(lefts) == pytest.approx(N / 2 * (1 + np.exp(-2 * k * t)), abs=3.0)


@pytest.mark.parametrize("bimolecular_side", ["forward", "reverse"])
def test_bimolecular_binding_matches_ode(bimolecular_side):
    """A + B <=> C in 1 fL: the Gillespie equilibrium mean of C matches the ODE (about 47 of 100)."""
    A, B, C = Species("A"), Species("B"), Species("C")
    kon, koff = 1e6 / (unit.molar * unit.sec), 0.1 / unit.sec
    if bimolecular_side == "forward":
        rxn = Reaction("bind", [A, B], [C], [1, 1], [1], kf=kon, kr=koff)
    else:
        rxn = Reaction("unbind", [C], [A, B], [1], [1, 1], kf=koff, kr=kon)
    flat = _flat([rxn], volume=1e-15 * unit.L)

    ode = ODESystem(flat)
    idx = ode.state.index["main"]
    ode.set_q([idx["A"], idx["B"]], 100)
    ode.run(60.0)
    c_ode = ode.state.q_val[idx["C"]]

    finals = []
    for _ in range(30):
        g = GillespieSystem(flat)
        gi = g.state.index["main"]
        g.set_q([gi["A"], gi["B"]], 100)
        g.run(30.0)
        finals.append(g.state.q_val[gi["C"]])
    assert np.mean(finals) == pytest.approx(c_ode, abs=3.0)


def test_decay_to_extinction_completes(decay_model):
    """Pure decay of 5 molecules over 50 half-lives ends at zero without error."""
    g = GillespieSystem(decay_model["flat"])
    i = g.state.index["main"]["A"]
    g.set_q([i], 5)
    g.run(70.0)
    assert g.state.q_val[i] == 0


def test_single_molecule_survival():
    """One A with A -> B at 1/s survives 0.1 s with probability ~exp(-0.1) = 0.90."""
    A, B = Species("A"), Species("B")
    kf, kr, t = 1.0, 1e-3, 0.1
    flat = _flat([Reaction("AB", [A], [B], [1], [1], kf=kf / unit.sec, kr=kr / unit.sec)])
    survived = []
    for _ in range(200):
        g = GillespieSystem(flat)
        i = g.state.index["main"]["A"]
        g.set_q([i], 1)
        g.run(t)
        survived.append(g.state.q_val[i])
    expected = (kr + kf * np.exp(-(kf + kr) * t)) / (kf + kr)
    assert np.mean(survived) == pytest.approx(expected, abs=0.07)
