"""Tests for volume-scaled transport in ODESystem: a compartment exchanging
with a Reservoir through a ResConnection (constant and time-dependent
reservoir concentrations), and DivByVConnection diffusion between two
compartments of unequal volume, each against its analytic solution.

Geometry used for the reservoir tests: a 2 x 1 x 1 nm box (V = 2 nm^3) with
a ResConnection on its x face, D = 1 nm^2/s, so kV = D * 1 nm^2 / 1 nm (half
the 2 nm width) = 1 nm^3/s and the exchange rate is k = kV / V = 0.5 /s. Reservoir
concentrations are in molecules/nm^3; the state holds molecule counts.
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.compartments.compartment import Compartment, Compartment3D, Reservoir
from openrxn.connections import ResConnection, DivByVConnection
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem

NM = unit.nm
KV, V, K = 1.0, 2.0, 0.5


def _reservoir_system(species, **res_kwargs):
    """ODESystem for the 2 nm^3 box 'c' coupled to Reservoir 'res' for each species."""
    c = Compartment3D("c", pos=[(0 * NM, 2 * NM), (0 * NM, 1 * NM), (0 * NM, 1 * NM)],
                      surface_area={"yz": 1 * NM ** 2, "xz": 2 * NM ** 2, "xy": 2 * NM ** 2})
    res = Reservoir("res", **res_kwargs)
    c.connect(res, ResConnection({s: 1.0 * NM ** 2 / unit.sec for s in species}, face="x"))
    return ODESystem(Model(compartments=[c, res]).flatten())


def test_constant_reservoir_fills_compartment():
    """With a fixed reservoir concentration C, n(t) = C V (1 - exp(-k t)) from an empty box."""
    C = 50.0
    sys = _reservoir_system(["A"], conc_funcs={"A": lambda t: C})
    assert sys.state.index["res"] == {}
    i = sys.state.index["c"]["A"]
    res = sys.run(20.0)
    assert res.y[i] == pytest.approx(C * V * (1 - np.exp(-K * res.t)), rel=1e-4, abs=1e-3)
    assert res.y[i][-1] == pytest.approx(C * V, rel=1e-2)


def test_time_dependent_reservoir_is_tracked_per_species():
    """Each species follows its own reservoir function: A sees a constant, B a decaying C0 exp(-a t)."""
    C_A, C0, a = 50.0, 100.0, 1.0
    sys = _reservoir_system(["A", "B"], conc_funcs={"A": lambda t: C_A, "B": lambda t: C0 * np.exp(-a * t)})
    iA, iB = sys.state.index["c"]["A"], sys.state.index["c"]["B"]
    sys.set_q([iA], 30.0)
    res = sys.run(15.0)
    t = res.t
    expected_A = C_A * V + (30.0 - C_A * V) * np.exp(-K * t)
    expected_B = KV * C0 / (K - a) * (np.exp(-a * t) - np.exp(-K * t))
    assert res.y[iA] == pytest.approx(expected_A, rel=1e-4, abs=1e-3)
    assert res.y[iB] == pytest.approx(expected_B, rel=1e-4, abs=1e-3)


def test_reservoir_concs_quantity_matches_conc_func():
    """Reservoir(concs={'A': 50 /nm^3}) gives the same filling curve as conc_funcs returning 50."""
    sys = _reservoir_system(["A"], concs={"A": 50.0 / NM ** 3})
    i = sys.state.index["c"]["A"]
    res = sys.run(20.0)
    assert res.y[i] == pytest.approx(50.0 * V * (1 - np.exp(-K * res.t)), rel=1e-4, abs=1e-3)


def test_divbyv_diffusion_between_unequal_volumes():
    """kV = 2/3 nm^3/s between 1 and 2 nm^3 boxes: n1 = 30 + 60 exp(-t), equal concentrations at the end."""
    left = Compartment("left", volume=1e-24 * unit.L)    # 1 nm^3
    right = Compartment("right", volume=2e-24 * unit.L)  # 2 nm^3
    conn = DivByVConnection({"A": 2.0 / 3.0 * NM ** 3 / unit.sec})
    left.connect(right, conn)
    right.connect(left, conn)
    sys = ODESystem(Model(compartments=[left, right]).flatten())
    iL, iR = sys.state.index["left"]["A"], sys.state.index["right"]["A"]
    sys.set_q([iL], 90.0)
    res = sys.run(10.0)
    # relaxation rate kV (1/V1 + 1/V2) = 1 /s; equilibrium n1/V1 = n2/V2 with n1 + n2 = 90
    assert res.y[iL] == pytest.approx(30.0 + 60.0 * np.exp(-res.t), rel=1e-4, abs=1e-3)
    assert res.y[iL] + res.y[iR] == pytest.approx(np.full(res.t.size, 90.0), rel=1e-6)
    assert res.y[iR][-1] / res.y[iL][-1] == pytest.approx(2.0, rel=1e-3)
