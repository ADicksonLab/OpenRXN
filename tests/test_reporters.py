"""Tests for openrxn.systems.reporters. Reporters are pure functions of an
input state vector, so they're exercised directly (no system run).
"""
import numpy as np
import pytest

from openrxn.systems.reporters import (
    Reporter, AllReporter, SelectionReporter, SumReporter, AvgReporter, MaxReporter, MinReporter,
)

VEC = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0])


def test_all_reporter():
    """AllReporter stores the full state vector with its timestamp."""
    r = AllReporter(freq=1)
    r.report(0.0, VEC)
    rep = r.reports()[0]
    assert rep["t"] == 0.0
    assert np.array_equal(rep["report"], VEC)


def test_selection_reporter():
    """SelectionReporter returns just the selected indices of the state vector."""
    r = SelectionReporter([1, 2, 4])
    r.report(0.0, VEC)
    assert np.array_equal(r.reports()[0]["report"], np.array([1.0, 4.0, 5.0]))


def test_sum_reporter():
    """SumReporter returns the sum over the selected indices."""
    r = SumReporter([1, 2, 4])
    r.report(0.0, VEC)
    assert r.reports()[0]["report"] == pytest.approx(10.0)


def test_avg_reporter():
    """AvgReporter returns the mean over the selected indices."""
    r = AvgReporter([1, 2, 4])
    r.report(0.0, VEC)
    assert r.reports()[0]["report"] == pytest.approx(10.0 / 3.0)


def test_max_reporter():
    """MaxReporter returns (max, argmax) within the selection."""
    r = MaxReporter([0, 1, 2])  # values [3, 1, 4]
    r.report(0.0, VEC)
    val, arg = r.reports()[0]["report"]
    assert val == pytest.approx(4.0)
    assert arg == 2


def test_min_reporter_returns_min():
    """MinReporter should return (min, argmin) within the selection."""
    r = MinReporter([0, 1, 2])  # values [3, 1, 4]
    r.report(0.0, VEC)
    val, arg = r.reports()[0]["report"]
    assert val == pytest.approx(1.0)
    assert arg == 1

def test_reports_accumulate_in_call_order():
    """Each report() call appends one {'t', 'report'} entry, in call order."""
    r = SumReporter([0, 5])
    r.report(0.5, VEC)
    r.report(1.5, 2 * VEC)
    assert [x["t"] for x in r.reports()] == [0.5, 1.5]
    assert [x["report"] for x in r.reports()] == pytest.approx([12.0, 24.0])


def test_reporter_frequencies():
    """Concrete reporters default to freq=1; the base Reporter defaults to None."""
    assert AllReporter().freq == 1
    assert SelectionReporter([0], freq=2.5).freq == 2.5
    assert Reporter().freq is None
    assert Reporter().reports() == []


def test_base_reporter_report_raises():
    """The abstract Reporter.report must raise NotImplementedError."""
    with pytest.raises(NotImplementedError):
        Reporter(freq=1).report(0.0, VEC)
