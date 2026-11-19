"""Tests for the fixes to the open OpenRXN issues (2026-11-04):

1. set_q accepts any amount (mol, mmol) or concentration (mol/L, molar, nM)
   and converts nm**3 volumes correctly.
2. ResConnection uses half the cell width (center to the shared face), which
   gives the exact steady-state profile between two reservoirs.
3. run() reports t = 0 and returns the whole trajectory when reporters split
   the run into segments.
4. GillespieSystem supports constant-concentration reservoirs.
5. Links between arrays use the reversed connection for the way back, so
   direction-dependent rates conserve mass; asymmetric links inside one array
   are rejected; Model.flatten warns about links that do not mirror.
"""
import logging

import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment, Compartment1D, Compartment3D, Reservoir
from openrxn.compartments.arrays import CompartmentArray1D, CompartmentArray3D
from openrxn.connections import IsotropicConnection, AnisotropicConnection, FicksConnection, ResConnection
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem
from openrxn.systems.GillespieSystem import GillespieSystem
from openrxn.systems.reporters import AllReporter

NM = unit.nm
NA = 6.022e23
PER_S = 1 / unit.sec


# --- 1. set_q units --------------------------------------------------------------

def _decay_sys(comp):
    A = Species("A")
    comp.add_rxn_to_compartment(Reaction("decay", [A], [], [1], [], kf=0.5 / unit.sec))
    sys = ODESystem(Model(compartments=[comp]).flatten())
    return sys, sys.state.index[comp.ID]["A"]


@pytest.mark.parametrize("Q, expected", [
    (1e-6 * unit.molar, 1e-6 * 2.0 * NA),
    (1.0 * unit.millimolar, 1e-3 * 2.0 * NA),
    (5.0 * unit.nanomolar, 5e-9 * 2.0 * NA),
    (1e-6 * unit.mol / unit.L, 1e-6 * 2.0 * NA),
    (2.0 * unit.mmol, 2e-3 * NA),
])
def test_set_q_accepts_any_amount_or_concentration_unit(Q, expected):
    """molar, mM, nM, mol/L (x 2 L volume x NA) and mmol (x NA) all give molecule counts."""
    sys, i = _decay_sys(Compartment("main", volume=2.0 * unit.L))
    sys.set_q([i], Q)
    assert sys.state.q_val[i] == pytest.approx(expected, rel=1e-12)


def test_set_q_mol_is_exactly_n_times_avogadro():
    """The mol path the paper notebooks use is unchanged: q = n * NA exactly."""
    sys, i = _decay_sys(Compartment("main", volume=1.0 * unit.L))
    sys.set_q([i], 3.7e-9 * unit.mol)
    assert sys.state.q_val[i] == 3.7e-9 * sys.NA


def test_set_q_concentration_in_nm3_volume():
    """1 mM in a 1000 nm^3 box is 1e-3 mol/L * 1e-21 L * NA = 0.6022 molecules (was 1e24 too many)."""
    sys, i = _decay_sys(Compartment("box", volume=1000 * NM ** 3))
    sys.set_q([i], 1.0 * unit.millimolar)
    assert sys.state.q_val[i] == pytest.approx(1e-3 * 1e-21 * NA, rel=1e-12)


def test_set_q_concentration_needs_3d_volume():
    """A 1D compartment's 'volume' is a length, so a concentration cannot be converted."""
    sys, i = _decay_sys(Compartment1D("seg", pos=[(0 * NM, 2 * NM)]))
    with pytest.raises(ValueError, match="3D volume"):
        sys.set_q([i], 1.0 * unit.molar)


def test_set_q_accepts_single_int_index():
    sys, i = _decay_sys(Compartment("main", volume=1.0 * unit.L))
    sys.set_q(i, 4.0)
    assert sys.state.q_val[i] == 4.0


# --- 2. ResConnection distance: half the cell width -------------------------------

def test_reservoir_half_cell_gives_exact_steady_profile():
    """N cells between reservoirs at C0 (x = 0) and 0 (x = L): at steady state each cell
    center holds C0 (1 - x/L), the exact continuum profile (a full-width distance would
    lower every value by a factor N/(N+1))."""
    N, L, C0, D = 5, 10.0, 2.0, 1.0
    arr = CompartmentArray3D("row", np.linspace(0, L, N + 1) * NM, np.array([0, 1.0]) * NM,
                             np.array([0, 1.0]) * NM, FicksConnection({"A": D * NM ** 2 / unit.sec}))
    left = Reservoir("left", concs={"A": C0 / NM ** 3})
    right = Reservoir("right", concs={"A": 0.0 / NM ** 3})
    arr.compartments[(0, 0, 0)].connect(left, ResConnection({"A": D * NM ** 2 / unit.sec}, face="x"))
    arr.compartments[(N - 1, 0, 0)].connect(right, ResConnection({"A": D * NM ** 2 / unit.sec}, face="x"))
    sys = ODESystem(Model(arrays=[arr], compartments=[left, right]).flatten())
    sys.run(600.0, rtol=1e-10, atol=1e-12)
    h = L / N
    for k in range(N):
        x = (k + 0.5) * h
        n = sys.state.q_val[sys.state.index[f"row-{k}_0_0"]["A"]]
        assert n / (h * 1.0 * 1.0) == pytest.approx(C0 * (1 - x / L), rel=1e-6)


# --- 3. run(): t = 0 report and the whole trajectory -----------------------------

def test_run_with_reporters_returns_whole_trajectory(decay_model):
    """With a freq = 1 reporter, run(5) returns t from 0 to 5 (no repeated points) and the
    analytic decay; the first report is the initial state at t = 0."""
    sys = ODESystem(decay_model["flat"])
    i = sys.state.index["main"]["A"]
    k = decay_model["kf"].magnitude
    sys.set_q([i], 100.0)
    rep = AllReporter(freq=1.0)
    sys.add_reporter(rep)
    res = sys.run(5.0)
    assert res.success
    assert res.t[0] == 0.0 and res.t[-1] == pytest.approx(5.0)
    assert np.all(np.diff(res.t) > 0)
    assert res.y[i] == pytest.approx(100.0 * np.exp(-k * res.t), rel=1e-4)
    assert rep.reports()[0]["t"] == 0
    assert rep.reports()[0]["report"][i] == 100.0


def test_reports_are_copies_of_the_state(decay_model):
    """Changing the state after run() does not change what the reporter recorded."""
    sys = ODESystem(decay_model["flat"])
    i = sys.state.index["main"]["A"]
    sys.set_q([i], 100.0)
    rep = AllReporter(freq=1.0)
    sys.add_reporter(rep)
    sys.run(2.0)
    before = rep.reports()[-1]["report"][i]
    sys.set_q([i], 7.0)
    assert rep.reports()[-1]["report"][i] == before


def test_run_with_reporters_and_t_eval(decay_model):
    """One t_eval passed to run() is split over the segments and joined back without duplicates."""
    sys = ODESystem(decay_model["flat"])
    i = sys.state.index["main"]["A"]
    sys.set_q([i], 100.0)
    sys.add_reporter(AllReporter(freq=1.0))
    t_eval = np.linspace(0, 5, 11)
    res = sys.run(5.0, t_eval=t_eval)
    assert np.array_equal(res.t, t_eval)
    assert res.y.shape == (sys.state.size, 11)


# --- 4. Gillespie with reservoirs -------------------------------------------------

def _res_box(system_cls, res_kwargs, conn=None):
    """2 x 1 x 1 nm box linked on its x face to a reservoir: kV = 1 nm^3/s, k = kV/V = 0.5 /s."""
    c = Compartment3D("c", pos=[(0 * NM, 2 * NM), (0 * NM, 1 * NM), (0 * NM, 1 * NM)],
                      surface_area={"yz": 1 * NM ** 2, "xz": 2 * NM ** 2, "xy": 2 * NM ** 2})
    res = Reservoir("res", **res_kwargs)
    c.connect(res, conn if conn is not None else ResConnection({"A": 1.0 * NM ** 2 / unit.sec}, face="x"))
    return system_cls(Model(compartments=[c, res]).flatten())


def test_gillespie_constant_reservoir_stationary_mean():
    """Box exchanging with a reservoir at C = 5 /nm^3: counts are Poisson with mean C V = 10."""
    g = _res_box(GillespieSystem, {"concs": {"A": 5.0 / NM ** 3}})
    i = g.state.index["c"]["A"]
    rep = AllReporter(freq=1.0)
    g.add_reporter(rep)
    g.run(4000.0)
    counts = np.array([r["report"][i] for r in rep.reports()][50:])
    assert counts.mean() == pytest.approx(10.0, abs=0.4)
    assert counts.var() == pytest.approx(10.0, rel=0.2)


def test_gillespie_reservoir_matches_ode_rates():
    """Gillespie uses the same in/out rates as the ODE: out k = 0.5 /s, in kV*C = 5 molecules/s."""
    g = _res_box(GillespieSystem, {"concs": {"A": 5.0 / NM ** 3}})
    i = g.state.index["c"]["A"]
    rates = sorted((p[0], len(p[1])) for p in g.processes)
    assert rates == [(pytest.approx(0.5), 1), (pytest.approx(5.0), 0)]
    assert all(d[0] == i for p in g.processes for d in p[2])


def test_gillespie_time_varying_reservoir_not_supported():
    with pytest.raises(NotImplementedError, match="time-varying"):
        _res_box(GillespieSystem, {"conc_funcs": {"A": lambda t: 5.0}})


def test_gillespie_reservoir_needs_res_connection():
    with pytest.raises(ValueError, match="ResConnection"):
        _res_box(GillespieSystem, {"concs": {"A": 5.0 / NM ** 3}},
                 conn=IsotropicConnection({"A": 1.0 * PER_S}))


def test_reservoir_copy_keeps_constant_concs():
    res = Reservoir("r", concs={"A": 3.0 / NM ** 3}, conc_funcs={"B": lambda t: t})
    assert res.copy().const_concs == {"A": 3.0}


# --- 5. Direction-dependent links -------------------------------------------------

def _pair_1d(link):
    a = CompartmentArray1D("a", np.array([0, 1.0]) * NM, IsotropicConnection({"D": 1.0 * PER_S}))
    b = CompartmentArray1D("b", np.array([1, 2.0]) * NM, IsotropicConnection({"D": 1.0 * PER_S}))
    a.stack(b, link)
    return a, b


def test_stack_anisotropic_conserves_mass_and_partitions():
    """a -> b at 1e-2 /s, b -> a at 1e-1 /s: the total stays 100 and n_a/n_b -> 10
    (it grew to 272 in 10 s when one object served both directions)."""
    a, b = _pair_1d(AnisotropicConnection({"D": (1e-2 * PER_S, 1e-1 * PER_S)}))
    sys = ODESystem(Model([a, b]).flatten())
    ia, ib = sys.state.index["a-0"]["D"], sys.state.index["b-0"]["D"]
    sys.set_q([ia], 100.0)
    res = sys.run(300.0, rtol=1e-10, atol=1e-10)
    assert res.y[[ia, ib]].sum(axis=0) == pytest.approx(100.0, rel=1e-8)
    assert res.y[ia, -1] == pytest.approx(100.0 * 10 / 11, rel=1e-6)


def test_stack_back_link_has_reversed_rates():
    a, b = _pair_1d(AnisotropicConnection({"D": (1e-2 * PER_S, 1e-1 * PER_S)}))
    k_ab = a.compartments[0].connections["b-0"][1].species_rates["D"]
    k_ba = b.compartments[0].connections["a-0"][1].species_rates["D"]
    assert (k_ba[0], k_ba[1]) == (k_ab[1], k_ab[0])


@pytest.mark.parametrize("side", ["x-", "x+", "y-", "y+", "z-", "z+"])
def test_join3d_back_link_has_reversed_rates(side):
    one = np.array([0, 1.0]) * NM
    a = CompartmentArray3D("a", one, one, one, IsotropicConnection({"D": 1.0 * PER_S}))
    b = CompartmentArray3D("b", one, one, one, IsotropicConnection({"D": 1.0 * PER_S}))
    a.join3D(b, AnisotropicConnection({"D": (1.0 * PER_S, 3.0 * PER_S)}), append_side=side)
    assert a.compartments[(0, 0, 0)].connections["b-0_0_0"][1].species_rates["D"][0] == 1.0 * PER_S
    assert b.compartments[(0, 0, 0)].connections["a-0_0_0"][1].species_rates["D"][0] == 3.0 * PER_S


def test_change_inter_connection_type_reverses_back_links():
    a, b = _pair_1d(IsotropicConnection({"D": 1.0 * PER_S}))
    a.change_all_inter_connection_type(b, AnisotropicConnection({"D": (2.0 * PER_S, 5.0 * PER_S)}))
    assert a.compartments[0].connections["b-0"][1].species_rates["D"][0] == 2.0 * PER_S
    assert b.compartments[0].connections["a-0"][1].species_rates["D"][0] == 5.0 * PER_S


@pytest.mark.parametrize("conn", [AnisotropicConnection({"D": (1.0 * PER_S, 2.0 * PER_S)}),
                                  IsotropicConnection({"D": (1.0 * PER_S, 2.0 * PER_S)})])
def test_asymmetric_link_inside_array_raises(conn):
    with pytest.raises(ValueError, match="inside an array"):
        CompartmentArray1D("a", np.arange(4.0) * NM, conn)


def test_symmetric_reverse_is_same_object():
    conn = IsotropicConnection({"D": 1.0 * PER_S})
    assert conn.reverse() is conn
    assert FicksConnection({"D": 1.0 * NM ** 2 / unit.sec}).reverse() is not None


def _two(c2_link):
    c1, c2 = Compartment("c1"), Compartment("c2")
    conn = AnisotropicConnection({"A": (1.0 * PER_S, 2.0 * PER_S)})
    c1.connect(c2, conn)
    if c2_link is not None:
        c2.connect(c1, c2_link(conn))
    return Model(compartments=[c1, c2])


@pytest.mark.parametrize("c2_link", [lambda conn: conn, None], ids=["same-object", "one-way"])
def test_flatten_warns_on_unmirrored_link(caplog, c2_link):
    with caplog.at_level(logging.WARNING):
        _two(c2_link).flatten()
    assert "do not mirror" in caplog.text


def test_flatten_quiet_on_mirrored_link(caplog):
    with caplog.at_level(logging.WARNING):
        flat = _two(lambda conn: conn.reverse()).flatten()
    assert "do not mirror" not in caplog.text
    assert Model.check_link_consistency(flat) == []


# --- 6. State is not shared with results; species order is deterministic --------

def test_set_q_after_run_keeps_previous_result(decay_model):
    """Pre-incubate, then change the state: the first result's last point stays as computed."""
    sys = ODESystem(decay_model["flat"])
    i = sys.state.index["main"]["A"]
    sys.set_q([i], 100.0)
    first = sys.run(1.0)
    last = first.y[i, -1]
    sys.set_q([i], 7.0)
    assert first.y[i, -1] == last
    sys.run(1.0)
    assert first.y[i, -1] == last


_ORDER_SCRIPT = """
from openrxn import unit
from openrxn.reactions import Reaction, Species
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.state import State
D, R, DR = Species('D'), Species('R'), Species('DR')
c = Compartment('reactor', volume=1 * unit.L)
c.add_rxn_to_compartment(Reaction('bind', [D, R], [DR], [1, 1], [1],
                                  kf=1e6 / (unit.molar * unit.sec), kr=0.1 / unit.sec))
print(','.join(State(model=Model(compartments=[c]).flatten()).species))
"""


def test_species_order_is_first_seen_and_same_in_every_process():
    """Species are indexed in first-seen order (reactants, then products), whatever PYTHONHASHSEED is."""
    import os, subprocess, sys as _sys
    orders = set()
    for seed in ["0", "1", "2", "3"]:
        env = dict(os.environ, PYTHONHASHSEED=seed)
        out = subprocess.run([_sys.executable, "-c", _ORDER_SCRIPT], env=env,
                             capture_output=True, text=True, check=True)
        orders.add(out.stdout.strip().splitlines()[-1])
    assert orders == {"D,R,DR"}
