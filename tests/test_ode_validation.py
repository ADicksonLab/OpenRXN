"""Validation tests for the ODE engine against analytic solutions:
first-order decay, A<=>B equilibrium, bimolecular mass conservation,
and two-compartment diffusion. Also set_q unit conversion and the volume
factor of bimolecular rates (kon / (NA V)) against the analytic Kd equilibrium.
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem


def test_first_order_decay(decay_model):
    """A --k--> : A(t) matches the closed form A0 * exp(-k t)."""
    sys = ODESystem(decay_model["flat"])
    i = sys.state.index["main"]["A"]
    sys.set_q([i], 100.0)
    res = sys.run(4.0)
    k = decay_model["kf"].to(1 / unit.sec).magnitude
    assert res.y[i][-1] == pytest.approx(100.0 * np.exp(-k * 4.0), rel=5e-3)


def test_equilibrium_ratio_and_conservation(equilibrium_model):
    """A<=>B reaches B/A = kf/kr while conserving total A + B."""
    sys = ODESystem(equilibrium_model["flat"])
    iA = sys.state.index["main"]["A"]
    iB = sys.state.index["main"]["B"]
    sys.set_q([iA], 100.0)
    sys.set_q([iB], 0.0)
    res = sys.run(50.0)
    A, B = res.y[iA][-1], res.y[iB][-1]
    kf = equilibrium_model["kf"].to(1 / unit.sec).magnitude
    kr = equilibrium_model["kr"].to(1 / unit.sec).magnitude
    assert A + B == pytest.approx(100.0, rel=1e-3)
    assert B / A == pytest.approx(kf / kr, rel=1e-2)


def test_bimolecular_mass_conservation():
    """A + B --> C conserves A + C (and C grows from zero)."""
    A, B, C = Species("A"), Species("B"), Species("C")
    r = Reaction("ABC", [A, B], [C], [1, 1], [1], kf=0.01 / unit.sec)
    c = Compartment("main")
    c.add_rxn_to_compartment(r)
    flat = Model(compartments=[c]).flatten()
    sys = ODESystem(flat)
    iA = sys.state.index["main"]["A"]
    iC = sys.state.index["main"]["C"]
    sys.set_q([iA], 50.0)
    sys.set_q([sys.state.index["main"]["B"]], 50.0)
    sys.set_q([iC], 0.0)
    res = sys.run(5.0)
    assert res.y[iA][-1] + res.y[iC][-1] == pytest.approx(50.0, rel=1e-3)
    assert res.y[iC][-1] > 0.0


def test_set_q_unit_conversion():
    """set_q turns mol and mol/L into molecule counts (x NA, x volume) and rejects other units."""
    A = Species("A")
    c = Compartment("main", volume=2.0 * unit.L)
    c.add_rxn_to_compartment(Reaction("decay", [A], [], [1], [], kf=0.5 / unit.sec))
    sys = ODESystem(Model(compartments=[c]).flatten())
    i = sys.state.index["main"]["A"]
    sys.set_q([i], 3e-9 * unit.mol)
    assert sys.state.q_val[i] == pytest.approx(3e-9 * sys.NA, rel=1e-12)
    sys.set_q([i], 3e-9 * unit.mol / unit.L)
    assert sys.state.q_val[i] == pytest.approx(3e-9 * 2.0 * sys.NA, rel=1e-12)
    sys.set_q([i], 7.0 * unit.dimensionless)
    assert sys.state.q_val[i] == 7.0
    with pytest.raises(ValueError):
        sys.set_q([i], 1.0 * unit.gram)


def test_bimolecular_equilibrium_with_volume():
    """D + R <=> DR in 2 L reaches the analytic Kd equilibrium; the second pair is written
    in reverse (DR <=> D + R) so the volume factor is also tested on kr."""
    kon, koff, V, D0, R0 = 1e8, 0.1, 2.0, 10e-9, 2e-9
    D1, R1, C1, D2, R2, C2 = (Species(n) for n in ["D1", "R1", "C1", "D2", "R2", "C2"])
    c = Compartment("main", volume=V * unit.L)
    c.add_rxns_to_compartment([
        Reaction("fwd", [D1, R1], [C1], [1, 1], [1], kf=kon / (unit.molar * unit.sec), kr=koff / unit.sec),
        Reaction("rev", [C2], [D2, R2], [1], [1, 1], kf=koff / unit.sec, kr=kon / (unit.molar * unit.sec))])
    sys = ODESystem(Model(compartments=[c]).flatten())
    idx = sys.state.index["main"]
    for s, conc in [("D1", D0), ("R1", R0), ("C1", 0.0), ("D2", D0), ("R2", R0), ("C2", 0.0)]:
        sys.set_q([idx[s]], conc * unit.mol / unit.L)
    res = sys.run(200.0)
    b = D0 + R0 + koff / kon
    C_eq = (b - np.sqrt(b ** 2 - 4 * D0 * R0)) / 2
    for s in ["C1", "C2"]:
        assert res.y[idx[s]][-1] / (sys.NA * V) == pytest.approx(C_eq, rel=1e-4)


def test_two_box_diffusion_equilibrates(two_box_model):
    """Symmetric diffusion between equal-volume boxes equilibrates to 50/50."""
    sys = ODESystem(two_box_model["flat"])
    iL = sys.state.index["left"]["A"]
    iR = sys.state.index["right"]["A"]
    sys.set_q([iL], 100.0)
    sys.set_q([iR], 0.0)
    res = sys.run(20.0)
    assert res.y[iL][-1] == pytest.approx(50.0, rel=1e-2)
    assert res.y[iR][-1] == pytest.approx(50.0, rel=1e-2)