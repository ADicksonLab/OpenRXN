"""Tests for openrxn.systems.system.System: the run() checkpoint loop, reporter
scheduling and values, kwargs forwarding, honouring a propagator's final_t,
init_state, and the abstract propagate. Also covers DerivFuncBuilder terms
(openrxn.systems.deriv) and the zero-order branch of the ODE Jacobian.

Deterministic runs use ODESystem on the conftest decay model, A(t) = A0 exp(-kf t).
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.model import Model
from openrxn.systems.system import System
from openrxn.systems.state import State
from openrxn.systems.ODESystem import ODESystem
from openrxn.systems.deriv import DerivFuncBuilder
from openrxn.systems.reporters import AllReporter, SelectionReporter, SumReporter

A0 = 100.0


class StubSystem(System):
    """System whose propagate adds 1 to every count, records its calls, and overshoots by self.overshoot."""

    overshoot = 0.25

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calls = []

    def propagate(self, t_interval, **kwargs):
        self.calls.append((t_interval, kwargs))
        self.state.q_val = self.state.q_val + 1
        return {"final_t": t_interval[1] + self.overshoot}


@pytest.fixture
def decay(decay_model):
    """ODESystem for A --0.5/s--> nothing with A0 = 100, plus the index of A and kf in 1/s."""
    sys = ODESystem(decay_model["flat"])
    i = sys.state.index["main"]["A"]
    sys.set_q([i], A0)
    return sys, i, decay_model["kf"].magnitude


def test_run_without_reporters_spans_whole_interval(decay):
    """With no reporters run() integrates [0, T] in one step and returns that solution."""
    sys, i, k = decay
    result = sys.run(4.0)
    assert result.t[0] == 0.0
    assert result.t[-1] == pytest.approx(4.0)
    assert sys.state.q_val[i] == pytest.approx(A0 * np.exp(-k * 4.0), rel=1e-4)


def test_add_reporter_reports_at_multiples_of_freq(decay):
    """An added AllReporter(freq=1) records the state at t = 0, 1, 2, 3 with the analytic values."""
    sys, i, k = decay
    rep = AllReporter(freq=1)
    sys.add_reporter(rep)
    sys.run(3.0)
    times = [r["t"] for r in rep.reports()]
    assert times == pytest.approx([0.0, 1.0, 2.0, 3.0])
    vals = [r["report"][i] for r in rep.reports()]
    assert vals == pytest.approx(A0 * np.exp(-k * np.array(times)), rel=1e-4)


def test_add_reporters_with_different_freqs(decay):
    """Reporters with freq 1 and 2.5 each report only at their own multiples, with correct values."""
    sys, i, k = decay
    every_1 = SumReporter([i], freq=1)
    every_2p5 = SelectionReporter([i], freq=2.5)
    sys.add_reporters([every_1, every_2p5])
    sys.run(5.0)
    t1 = [r["t"] for r in every_1.reports()]
    t2 = [r["t"] for r in every_2p5.reports()]
    assert t1 == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    assert t2 == pytest.approx([0.0, 2.5, 5.0])
    assert [r["report"] for r in every_1.reports()] == pytest.approx(A0 * np.exp(-k * np.array(t1)), rel=1e-4)
    assert [r["report"][0] for r in every_2p5.reports()] == pytest.approx(A0 * np.exp(-k * np.array(t2)), rel=1e-4)


def test_total_time_not_multiple_of_freq(decay):
    """With freq=2 and T=5 reports land at 0, 2 and 4, and the state still reaches t = 5."""
    sys, i, k = decay
    rep = AllReporter(freq=2)
    sys.add_reporter(rep)
    sys.run(5.0)
    assert [r["t"] for r in rep.reports()] == pytest.approx([0.0, 2.0, 4.0])
    assert sys.state.q_val[i] == pytest.approx(A0 * np.exp(-k * 5.0), rel=1e-4)


def test_constructor_reporters_are_copied(decay_model):
    """The reporters list given to System is copied, and systems never share a default list."""
    given = [AllReporter(freq=1)]
    s1 = System(decay_model["flat"], reporters=given)
    s1.add_reporter(AllReporter(freq=2))
    assert len(given) == 1
    assert len(s1.reporters) == 2
    s2 = System(decay_model["flat"])
    s2.add_reporter(AllReporter(freq=1))
    assert System(decay_model["flat"]).reporters == []


def test_run_forwards_kwargs_and_honours_final_t(decay_model):
    """run() passes kwargs to propagate and restarts each interval at the final_t it returns."""
    sys = StubSystem(decay_model["flat"])
    rep = AllReporter(freq=1)
    sys.add_reporter(rep)
    out = sys.run(3.0, method="X")
    assert [c[0] for c in sys.calls] == [(0, 1.0), (1.25, 2.0), (2.25, 3.0)]
    assert all(c[1] == {"method": "X"} for c in sys.calls)
    assert [r["t"] for r in rep.reports()] == pytest.approx([0.0, 1.25, 2.25, 3.25])
    i = sys.state.index["main"]["A"]
    assert [r["report"][i] for r in rep.reports()] == [0.0, 1.0, 2.0, 3.0]
    assert out == {"final_t": 3.25}


@pytest.mark.parametrize("total_time", [1 / 3, 0.1 + 0.2, 3e4 * 3**5, 3e4 * 3**11, 1e-7 / 3, 1800.0, 7])
def test_run_without_reporters_is_one_exact_call(decay_model, total_time):
    """Without reporters run() makes exactly one propagate call over (0, total_time), total_time untouched."""
    sys = StubSystem(decay_model["flat"])
    sys.run(total_time, rtol=1e-9)
    assert len(sys.calls) == 1
    (t0, t1), kwargs = sys.calls[0]
    assert t0 == 0 and type(t0) is int
    assert t1 == total_time and type(t1) is type(total_time)
    assert kwargs == {"rtol": 1e-9}


@pytest.mark.parametrize("total_time, freqs, n_reports", [(0.3, [0.1], [3]),
                                                          (1.0, [0.1], [10]),
                                                          (0.7, [0.1, 0.25], [7, 2])])
def test_run_with_reporters_keeps_endpoints_exact(decay_model, total_time, freqs, n_reports):
    """Round-off duplicates are merged: intervals chain from exactly 0 to exactly total_time, one report per multiple (and t = 0)."""
    sys = StubSystem(decay_model["flat"])
    sys.overshoot = 0.0
    reps = [AllReporter(freq=f) for f in freqs]
    sys.add_reporters(reps)
    sys.run(total_time)
    ivals = [c[0] for c in sys.calls]
    assert ivals[0][0] == 0
    assert ivals[-1][1] == total_time
    assert all(a[1] == b[0] for a, b in zip(ivals, ivals[1:]))
    assert all(t1 - t0 > 1e-6 for t0, t1 in ivals)
    for rep, f, n in zip(reps, freqs, n_reports):
        assert [r["t"] for r in rep.reports()] == pytest.approx([f * k for k in range(0, n + 1)])


def test_init_state_is_used(decay_model):
    """A State passed as init_state becomes system.state and is the starting point of run()."""
    flat = decay_model["flat"]
    st = State(model=flat)
    i = st.index["main"]["A"]
    st.q_val[i] = 40.0
    sys = ODESystem(flat, init_state=st)
    assert sys.state is st
    sys.run(2.0)
    assert sys.state.q_val[i] == pytest.approx(40.0 * np.exp(-decay_model["kf"].magnitude * 2.0), rel=1e-4)


def test_base_system_builds_state_and_propagate_is_abstract(decay_model):
    """The base System builds a State from the model; its propagate raises NotImplementedError."""
    sys = System(decay_model["flat"])
    assert sys.state.index == {"main": {"A": 0}}
    assert sys.reporters == []
    with pytest.raises(NotImplementedError):
        sys.propagate((0, 1.0))


def test_base_run_raises_not_implemented(decay_model):
    """Running a bare System reports the missing propagate as NotImplementedError."""
    with pytest.raises(NotImplementedError):
        System(decay_model["flat"]).run(1.0)


def test_reporter_survives_float_roundoff(decay):
    """freq=0.1 over T=0.3 reports at 0, 0.1, 0.2 and 0.3."""
    sys, i, k = decay
    rep = AllReporter(freq=0.1)
    sys.add_reporter(rep)
    sys.run(0.3)
    assert [r["t"] for r in rep.reports()] == pytest.approx([0.0, 0.1, 0.2, 0.3])


def test_zero_order_synthesis_ode_and_jacobian():
    """0 -> A at 2/s with A -> 0 at 0.5/s gives A(t) = 4 (1 - exp(-0.5 t)); the source adds nothing to J."""
    A = Species("A")
    c = Compartment("main")
    c.add_rxns_to_compartment([Reaction("syn", [], [A], [], [1], kf=2.0 / unit.sec),
                               Reaction("deg", [A], [], [1], [], kf=0.5 / unit.sec)])
    sys = ODESystem(Model(compartments=[c]).flatten())
    i = sys.state.index["main"]["A"]
    sys.run(3.0)
    assert sys.state.q_val[i] == pytest.approx(4.0 * (1 - np.exp(-1.5)), rel=1e-4)
    assert sys._jac(0.0, np.array([10.0])) == pytest.approx(np.array([[-0.5]]))
    assert sys._dQ_dt(0.0, np.array([10.0])) == pytest.approx([2.0 - 5.0])


def test_deriv_func_builder_terms():
    """DerivFuncBuilder converts rates to 1/s, scales by n_per, drops zero rates, and adds reservoir inflow."""
    sources = [(6.0 / unit.minute, [0], 1), (1.0 / unit.sec, [0, 1], 2)]
    sinks = [(0.5 / unit.sec, [1], 1), (0.0 / unit.sec, [0], 1)]
    reservoir = [(3.0 * unit.nm**3 / unit.sec, lambda t: 2.0 + t, 2)]
    d = DerivFuncBuilder(sources, sinks, reservoir)
    assert [(p, list(ix)) for p, ix in d.terms] == [(pytest.approx(0.1), [0]), (2.0, [0, 1]), (-0.5, [1])]
    Q = np.array([4.0, 5.0])
    assert d.deriv_func(Q, 1.0) == pytest.approx(0.1 * 4 + 2 * 4 * 5 - 0.5 * 5 + 6 * (2.0 + 1.0))
