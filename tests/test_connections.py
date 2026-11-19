"""Tests for openrxn.connections: rate normalization, directionality,
reverse(), __repr__, and Fick's-law resolution of FicksConnection and
ResConnection into a DivByVConnection.
"""
import logging

import pytest

from openrxn import unit
from openrxn.connections import (IsotropicConnection, AnisotropicConnection, DivByVConnection,
                                 FicksConnection, ResConnection)


def test_isotropic_symmetric():
    """A scalar isotropic rate is stored as a symmetric (k, k) tuple in 1/s."""
    conn = IsotropicConnection({"A": 0.5 / unit.sec})
    k_out, k_in = conn.species_rates["A"]
    assert k_out.to(1 / unit.sec).magnitude == pytest.approx(0.5)
    assert k_in.to(1 / unit.sec).magnitude == pytest.approx(0.5)


def test_anisotropic_keeps_direction():
    """An anisotropic (k_out, k_in) tuple keeps its two distinct rates."""
    conn = AnisotropicConnection({"A": (0.5 / unit.sec, 0.2 / unit.sec)})
    k_out, k_in = conn.species_rates["A"]
    assert k_out.to(1 / unit.sec).magnitude == pytest.approx(0.5)
    assert k_in.to(1 / unit.sec).magnitude == pytest.approx(0.2)


def test_anisotropic_reverse():
    """reverse() swaps the out/in rates of an anisotropic connection."""
    conn = AnisotropicConnection({"A": (0.5 / unit.sec, 0.2 / unit.sec)})
    rev = conn.reverse()
    k_out, k_in = rev.species_rates["A"]
    assert k_out.to(1 / unit.sec).magnitude == pytest.approx(0.2)
    assert k_in.to(1 / unit.sec).magnitude == pytest.approx(0.5)


def test_anisotropic_bad_tuple_length():
    """A rate tuple that is not length 2 raises ValueError."""
    with pytest.raises(ValueError):
        AnisotropicConnection({"A": (0.5 / unit.sec, 0.2 / unit.sec, 0.1 / unit.sec)})


def test_anisotropic_scalar_becomes_symmetric_with_warning(caplog):
    """A scalar anisotropic rate warns and is stored as (k, k), converted to 1/s."""
    with caplog.at_level(logging.WARNING):
        conn = AnisotropicConnection({"A": 6.0 / unit.minute})
    k_out, k_in = conn.species_rates["A"]
    assert k_out.units == k_in.units == unit.sec ** -1
    assert (k_out.magnitude, k_in.magnitude) == pytest.approx((0.1, 0.1))
    assert "Species A: one scalar rate provided" in caplog.text


def test_anisotropic_repr():
    """repr lists dim, species count, (k_out, k_in) magnitudes and the 1/s unit."""
    text = repr(AnisotropicConnection({"A": (0.5 / unit.sec, 0.2 / unit.sec)}, dim=2))
    assert text.startswith("AnisotropicConnection(dim=2,n_species=1,")
    assert "'A': (0.5, 0.2)" in text
    assert text.endswith("1 / second)")


def test_isotropic_repr_converts_units():
    """repr of an isotropic connection shows every species' rate after conversion to 1/s."""
    text = repr(IsotropicConnection({"A": 0.5 / unit.sec, "B": 6.0 / unit.minute}))
    assert text.startswith("IsotropicConnection(dim=3,n_species=2,")
    assert "'A': (0.5, 0.5)" in text
    assert "'B': (0.1, 0.1)" in text
    assert text.endswith("1 / second)")


@pytest.mark.parametrize("cls", [AnisotropicConnection, IsotropicConnection])
def test_repr_of_empty_connection(cls):
    """repr of a connection with no species returns a string instead of raising."""
    assert "n_species=0" in repr(cls({}))


def test_ficks_resolve_rate_and_units():
    """FicksConnection.resolve() yields a DivByVConnection with rate D*A/dx in nm^3/s."""
    D = 2.0 * unit.nm ** 2 / unit.sec
    conn = FicksConnection({"A": D}, surface_area=3.0 * unit.nm ** 2, ic_distance=6.0 * unit.nm)
    resolved = conn.resolve()
    assert isinstance(resolved, DivByVConnection)
    k_out, _ = resolved.species_rates["A"]
    # D*A/dx = 2 * 3 / 6 = 1 nm^3/s
    assert k_out.to(unit.nm ** 3 / unit.sec).magnitude == pytest.approx(1.0)


def test_ficks_resolve_unready_raises():
    """Resolving a FicksConnection without surface_area/ic_distance raises ValueError."""
    conn = FicksConnection({"A": 2.0 * unit.nm ** 2 / unit.sec})
    with pytest.raises(ValueError):
        conn.resolve()


def test_res_connection_stores_geometry():
    """ResConnection keeps the D constants, interface geometry, dim and face it was given."""
    D = {"A": 1.0 * unit.nm ** 2 / unit.sec}
    conn = ResConnection(D, surface_area=2.0 * unit.nm ** 2, ic_distance=4.0 * unit.nm, dim=3, face="z")
    assert conn.species_d_constants is D
    assert conn.surface_area == 2.0 * unit.nm ** 2
    assert conn.ic_distance == 4.0 * unit.nm
    assert (conn.dim, conn.face) == (3, "z")


def test_res_connection_resolve_rate_and_units():
    """ResConnection.resolve() gives a symmetric DivByVConnection with D*A/dx in nm^3/s, per species."""
    conn = ResConnection({"A": 1e-14 * unit.cm ** 2 / unit.sec, "B": 4.0 * unit.nm ** 2 / unit.sec},
                         surface_area=3.0 * unit.nm ** 2, ic_distance=6.0 * unit.nm)
    resolved = conn.resolve()
    assert isinstance(resolved, DivByVConnection)
    assert resolved.dim == 3
    # A: 1 nm^2/s * 3 / 6 = 0.5;  B: 4 * 3 / 6 = 2
    for s, kv in [("A", 0.5), ("B", 2.0)]:
        k_out, k_in = resolved.species_rates[s]
        assert k_out.units == unit.nm ** 3 / unit.sec
        assert (k_out.magnitude, k_in.magnitude) == pytest.approx((kv, kv))


def test_res_connection_resolve_2d_units():
    """In 2D, D*A/dx is expressed in nm^2/s (A is an interface length)."""
    conn = ResConnection({"A": 2.0 * unit.nm ** 2 / unit.sec}, surface_area=3.0 * unit.nm,
                         ic_distance=1.5 * unit.nm, dim=2)
    k_out, _ = conn.resolve().species_rates["A"]
    assert k_out.units == unit.nm ** 2 / unit.sec
    assert k_out.magnitude == pytest.approx(4.0)


@pytest.mark.parametrize("missing", ["surface_area", "ic_distance"])
def test_res_connection_resolve_unready_raises(missing):
    """Resolving a ResConnection without surface_area or ic_distance raises ValueError."""
    geom = {"surface_area": 1.0 * unit.nm ** 2, "ic_distance": 1.0 * unit.nm}
    geom[missing] = None
    conn = ResConnection({"A": 1.0 * unit.nm ** 2 / unit.sec}, **geom)
    with pytest.raises(ValueError, match="not ready"):
        conn.resolve()
