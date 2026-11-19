"""Tests for openrxn.compartments.arrays: grid geometry (positions, volumes,
face areas, box lengths), nearest-neighbor topology with and without periodic
wrapping, stack()/join3D() between arrays, array-wide reaction and
connection-type updates, to_graph export, and transport on a periodic ring against the
analytic solution.
"""
import networkx as nx
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.arrays import CompartmentArray1D, CompartmentArray2D, CompartmentArray3D
from openrxn.connections import IsotropicConnection, FicksConnection
from openrxn.model import Model
from openrxn.systems.ODESystem import ODESystem


def _nm(*vals):
    return np.array(vals, dtype=float) * unit.nm


def _grid(n):
    return np.arange(n + 1, dtype=float) * unit.nm


def _iso(k=1.0):
    return IsotropicConnection({"A": k / unit.sec})


def _tag(array_ID, idx):
    return array_ID + "-" + "_".join(str(i) for i in idx)


def _n_links(arr):
    return sum(len(c.connections) for c in arr.compartments.values())


@pytest.mark.parametrize("build", [
    lambda conn: CompartmentArray1D("a", _grid(2), conn),
    lambda conn: CompartmentArray2D("a", _grid(2), _grid(2), conn),
    lambda conn: CompartmentArray3D("a", _grid(2), _grid(2), _grid(2), conn),
], ids=["1D", "2D", "3D"])
def test_constructor_rejects_non_connection(build):
    """Every array constructor asserts conn_type is a Connection instance."""
    with pytest.raises(AssertionError):
        build({"A": 1.0 / unit.sec})


def test_1d_geometry():
    """1D cell i spans [p_i, p_i+1]; its volume is the width; box_len is the full span."""
    arr = CompartmentArray1D("m", _nm(0, 1, 3, 6), _iso())
    assert arr.n_compartments == 3
    assert sorted(arr.compartments) == [0, 1, 2]
    widths = [arr.compartments[i].volume.to(unit.nm).magnitude for i in range(3)]
    assert widths == pytest.approx([1.0, 2.0, 3.0])
    lo, hi = arr.compartments[1].pos[0]
    assert (lo.to(unit.nm).magnitude, hi.to(unit.nm).magnitude) == (1.0, 3.0)
    assert arr.box_len[0].to(unit.nm).magnitude == pytest.approx(6.0)
    assert all(c.array_ID == "m" for c in arr.compartments.values())


def test_1d_chain_neighbors():
    """A non-periodic 1D array links each cell to its immediate neighbors only."""
    conn = _iso()
    arr = CompartmentArray1D("m", _grid(3), conn)
    assert set(arr.compartments[0].connections) == {"m-1"}
    assert set(arr.compartments[1].connections) == {"m-0", "m-2"}
    assert set(arr.compartments[2].connections) == {"m-1"}
    assert arr.compartments[1].connections["m-2"] == (arr.compartments[2], conn)


@pytest.mark.parametrize("periodic, stored, n_links", [
    (False, [False], 6), (True, [True], 8), ([False], [False], 6), ([True], [True], 8),
])
def test_1d_periodic_stored_as_list(periodic, stored, n_links):
    """1D periodic is stored as a one-element list; bool and list inputs wrap the same way."""
    arr = CompartmentArray1D("m", _grid(4), _iso(), periodic=periodic)
    assert arr.periodic == stored
    assert _n_links(arr) == n_links


def test_1d_periodic_ring_matches_analytic():
    """Transport on a periodic 4-ring from a point source matches the circulant solution."""
    k, t = 1.0, 0.5
    arr = CompartmentArray1D("m", _grid(4), _iso(k), periodic=True)
    sys = ODESystem(Model(arrays=[arr]).flatten())
    idx = [sys.state.index[f"m-{i}"]["A"] for i in range(4)]
    sys.set_q([idx[0]], 100.0)
    sys.set_q(idx[1:], 0.0)
    res = sys.run(t)
    e2, e4 = np.exp(-2 * k * t), np.exp(-4 * k * t)
    expected = [25 * (1 + 2 * e2 + e4), 25 * (1 - e4), 25 * (1 - 2 * e2 + e4), 25 * (1 - e4)]
    assert [res.y[i][-1] for i in idx] == pytest.approx(expected, rel=5e-3)


def test_1d_stack_links_matching_indices():
    """stack() links cell i of each 1D array to cell i of the other, in both directions."""
    a = CompartmentArray1D("a", _grid(3), _iso())
    b = CompartmentArray1D("b", _grid(3), _iso())
    inter = _iso(0.2)
    a.stack(b, inter)
    for i in range(3):
        assert a.compartments[i].connections[f"b-{i}"] == (b.compartments[i], inter)
        assert b.compartments[i].connections[f"a-{i}"] == (a.compartments[i], inter)
    assert set(a.compartments[1].connections) == {"a-0", "a-2", "b-1"}


@pytest.mark.parametrize("make_other, conn", [
    (lambda: CompartmentArray2D("b", _grid(3), _grid(1), _iso()), _iso()),
    (lambda: CompartmentArray1D("b", _grid(2), _iso()), _iso()),
    (lambda: CompartmentArray1D("b", _grid(3), _iso()), "not a connection"),
], ids=["not_1d", "size_mismatch", "bad_conn"])
def test_1d_stack_validates_arguments(make_other, conn):
    """stack() rejects a non-1D array, a different length, or a non-Connection, linking nothing."""
    a = CompartmentArray1D("a", _grid(3), _iso())
    with pytest.raises(AssertionError):
        a.stack(make_other(), conn)
    assert _n_links(a) == 4


def test_2d_geometry():
    """2D cell (i,j) spans [x_i,x_i+1] x [y_j,y_j+1]; volume is its area; box_len is the extent."""
    arr = CompartmentArray2D("g", _nm(0, 1, 3), _nm(0, 2, 5, 6), _iso())
    assert (arr.nx, arr.ny, arr.n_compartments) == (2, 3, 6)
    assert set(arr.compartments) == {(i, j) for i in range(2) for j in range(3)}
    c = arr.compartments[(1, 1)]
    assert [(lo.magnitude, hi.magnitude) for lo, hi in c.pos] == [(1.0, 3.0), (2.0, 5.0)]
    assert c.volume.to(unit.nm ** 2).magnitude == pytest.approx(6.0)
    assert [L.to(unit.nm).magnitude for L in arr.box_len] == pytest.approx([3.0, 6.0])


def test_2d_grid_neighbors():
    """In a 3x3 grid a corner, edge and center cell have 2, 3 and 4 neighbors."""
    arr = CompartmentArray2D("g", _grid(3), _grid(3), _iso())
    assert set(arr.compartments[(0, 0)].connections) == {"g-1_0", "g-0_1"}
    assert set(arr.compartments[(0, 1)].connections) == {"g-1_1", "g-0_0", "g-0_2"}
    assert set(arr.compartments[(1, 1)].connections) == {"g-0_1", "g-2_1", "g-1_0", "g-1_2"}
    assert _n_links(arr) == 24


@pytest.mark.parametrize("periodic, wraps", [
    ([True, False], [(2, 0)]),
    ([False, True], [(0, 2)]),
    ([True, True], [(2, 0), (0, 2)]),
])
def test_2d_periodic_wraps(periodic, wraps):
    """Each periodic axis links the first and last cell of every row/column, both ways."""
    arr = CompartmentArray2D("g", _grid(3), _grid(3), _iso(), periodic=periodic)
    corner = arr.compartments[(0, 0)]
    assert set(corner.connections) == {"g-1_0", "g-0_1"} | {_tag("g", w) for w in wraps}
    for w in wraps:
        assert "g-0_0" in arr.compartments[w].connections
    assert _n_links(arr) == 24 + 6 * sum(periodic)


def test_2d_stack_links_matching_cells():
    """stack() links cell (i,j) of each 2D array to cell (i,j) of the other, in both directions."""
    a = CompartmentArray2D("a", _grid(2), _grid(3), _iso())
    b = CompartmentArray2D("b", _grid(2), _grid(3), _iso())
    inter = _iso(0.2)
    a.stack(b, inter)
    for idx, ca in a.compartments.items():
        cb = b.compartments[idx]
        assert ca.connections[_tag("b", idx)] == (cb, inter)
        assert cb.connections[_tag("a", idx)] == (ca, inter)
        assert sum(t.startswith("b-") for t in ca.connections) == 1


@pytest.mark.parametrize("make_other, conn", [
    (lambda: CompartmentArray1D("b", _grid(2), _iso()), _iso()),
    (lambda: CompartmentArray2D("b", _grid(3), _grid(2), _iso()), _iso()),
    (lambda: CompartmentArray2D("b", _grid(2), _grid(3), _iso()), _iso()),
    (lambda: CompartmentArray2D("b", _grid(2), _grid(2), _iso()), "not a connection"),
], ids=["not_2d", "x_mismatch", "y_mismatch", "bad_conn"])
def test_2d_stack_validates_arguments(make_other, conn):
    """stack() rejects a non-2D array, an x or y size mismatch, or a non-Connection."""
    a = CompartmentArray2D("a", _grid(2), _grid(2), _iso())
    with pytest.raises(AssertionError):
        a.stack(make_other(), conn)
    assert _n_links(a) == 8


def test_3d_geometry_and_face_areas():
    """3D cell volumes and face areas ('xy', 'yz', 'xz') follow from the per-axis widths."""
    arr = CompartmentArray3D("g", _nm(0, 1, 3), _nm(0, 2, 5), _nm(0, 4, 7), _iso())
    assert (arr.nx, arr.ny, arr.nz, arr.n_compartments) == (2, 2, 2, 8)
    assert [L.to(unit.nm).magnitude for L in arr.box_len] == pytest.approx([3.0, 5.0, 7.0])
    c = arr.compartments[(0, 1, 0)]
    assert [(lo.magnitude, hi.magnitude) for lo, hi in c.pos] == [(0.0, 1.0), (2.0, 5.0), (0.0, 4.0)]
    assert c.volume.to(unit.nm ** 3).magnitude == pytest.approx(12.0)
    areas = {f: a.to(unit.nm ** 2).magnitude for f, a in c.surface_area.items()}
    assert areas == pytest.approx({"xy": 3.0, "yz": 12.0, "xz": 4.0})
    c = arr.compartments[(1, 0, 1)]
    assert c.volume.to(unit.nm ** 3).magnitude == pytest.approx(12.0)
    areas = {f: a.to(unit.nm ** 2).magnitude for f, a in c.surface_area.items()}
    assert areas == pytest.approx({"xy": 4.0, "yz": 6.0, "xz": 6.0})


def test_3d_grid_neighbor_counts():
    """In a 3x3x3 grid each cell links to 3 (corner) up to 6 (center) face neighbors."""
    arr = CompartmentArray3D("g", _grid(3), _grid(3), _grid(3), _iso())
    for (i, j, k), c in arr.compartments.items():
        assert len(c.connections) == sum(2 if x == 1 else 1 for x in (i, j, k))
    assert set(arr.compartments[(1, 1, 1)].connections) == {
        "g-0_1_1", "g-2_1_1", "g-1_0_1", "g-1_2_1", "g-1_1_0", "g-1_1_2"}
    assert _n_links(arr) == 108


@pytest.mark.parametrize("periodic, wraps", [
    ([True, False, False], [(2, 0, 0)]),
    ([False, True, False], [(0, 2, 0)]),
    ([False, False, True], [(0, 0, 2)]),
    ([True, True, True], [(2, 0, 0), (0, 2, 0), (0, 0, 2)]),
])
def test_3d_periodic_wraps(periodic, wraps):
    """Each periodic axis links opposite boundary cells of the 3D grid, both ways."""
    arr = CompartmentArray3D("g", _grid(3), _grid(3), _grid(3), _iso(), periodic=periodic)
    corner = arr.compartments[(0, 0, 0)]
    assert set(corner.connections) == {"g-1_0_0", "g-0_1_0", "g-0_0_1"} | {_tag("g", w) for w in wraps}
    for w in wraps:
        assert "g-0_0_0" in arr.compartments[w].connections
    assert _n_links(arr) == 108 + 18 * sum(periodic)


@pytest.mark.parametrize("side", ["x-", "x+", "y-", "y+", "z-", "z+"])
def test_join3d_links_adjoining_faces(side):
    """join3D links self's `side` face to the opposite face of other, cell by cell, both ways."""
    axis = "xyz".index(side[0])
    other_shape = [2, 2, 2]
    other_shape[axis] = 3
    a = CompartmentArray3D("a", _grid(2), _grid(2), _grid(2), _iso())
    b = CompartmentArray3D("b", *[_grid(n) for n in other_shape], _iso())
    inter = _iso(0.2)
    a.join3D(b, inter, append_side=side)
    a_face = 0 if side[1] == "-" else 1
    b_face = 2 if side[1] == "-" else 0
    for idx, ca in a.compartments.items():
        cross = {t: v for t, v in ca.connections.items() if t.startswith("b-")}
        if idx[axis] != a_face:
            assert cross == {}
            continue
        b_idx = tuple(b_face if d == axis else idx[d] for d in range(3))
        cb = b.compartments[b_idx]
        assert cross == {_tag("b", b_idx): (cb, inter)}
        assert cb.connections[_tag("a", idx)] == (ca, inter)
    assert sum(t.startswith("a-") for c in b.compartments.values() for t in c.connections) == 4


@pytest.mark.parametrize("side", ["x+", "y-", "z+"])
def test_join3d_face_mismatch_raises(side):
    """join3D raises ValueError when the joined faces have different cell counts."""
    a = CompartmentArray3D("a", _grid(2), _grid(2), _grid(2), _iso())
    b = CompartmentArray3D("b", _grid(3), _grid(3), _grid(3), _iso())
    with pytest.raises(ValueError, match="dimensions of arrays don't match"):
        a.join3D(b, _iso(), append_side=side)
    assert _n_links(a) == 24


@pytest.mark.parametrize("make_other, conn, side", [
    (lambda: CompartmentArray2D("b", _grid(2), _grid(2), _iso()), _iso(), "x+"),
    (lambda: CompartmentArray3D("b", _grid(2), _grid(2), _grid(2), _iso()), _iso(), None),
    (lambda: CompartmentArray3D("b", _grid(2), _grid(2), _grid(2), _iso()), "not a connection", "x+"),
], ids=["not_3d", "missing_side", "bad_conn"])
def test_join3d_validates_arguments(make_other, conn, side):
    """join3D rejects a non-3D array, a missing/invalid append_side, or a non-Connection."""
    a = CompartmentArray3D("a", _grid(2), _grid(2), _grid(2), _iso())
    with pytest.raises(AssertionError):
        a.join3D(make_other(), conn, append_side=side)
    assert _n_links(a) == 24


def test_add_rxns_to_array_reaches_every_compartment():
    """add_rxns_to_array puts each reaction once in every compartment; repeats are ignored."""
    A, B = Species("A"), Species("B")
    r1 = Reaction("r1", [A], [B], [1], [1], kf=1.0 / unit.sec)
    r2 = Reaction("r2", [B], [], [1], [], kf=2.0 / unit.sec)
    arr = CompartmentArray2D("g", _grid(2), _grid(2), _iso())
    arr.add_rxns_to_array([r1, r2])
    arr.add_rxn_to_array(r1)
    for c in arr.compartments.values():
        assert [r.ID for r in c.reactions] == ["r1", "r2"]


def _stacked_pair():
    old, inter = _iso(1.0), _iso(0.2)
    a = CompartmentArray1D("a", _grid(3), old)
    b = CompartmentArray1D("b", _grid(3), old)
    a.stack(b, inter)
    return a, b, old, inter


def test_change_intra_connection_type_skips_inter_links():
    """change_all_intra_connection_type swaps only links whose neighbor is in the same array."""
    a, b, old, inter = _stacked_pair()
    new = _iso(5.0)
    a.change_all_intra_connection_type(new)
    for c in a.compartments.values():
        for tag, (_, conn) in c.connections.items():
            assert conn is (new if tag.startswith("a-") else inter)
    for c in b.compartments.values():
        for tag, (_, conn) in c.connections.items():
            assert conn is (old if tag.startswith("b-") else inter)
    assert _n_links(a) == 7


def test_change_inter_connection_type_skips_intra_links():
    """change_all_inter_connection_type swaps self->other links and keeps intra-array links."""
    a, b, old, inter = _stacked_pair()
    new = _iso(5.0)
    a.change_all_inter_connection_type(b, new)
    for c in a.compartments.values():
        for tag, (_, conn) in c.connections.items():
            assert conn is (new if tag.startswith("b-") else old)
    for c in b.compartments.values():
        for tag, (_, conn) in c.connections.items():
            assert conn is (new if tag.startswith("a-") else old)
    assert _n_links(a) == _n_links(b) == 7


def test_change_inter_connection_type_conserves_mass():
    """After changing the inter-array rate, A exchanged between two stacked cells is conserved."""
    a = CompartmentArray1D("a", _grid(1), _iso())
    b = CompartmentArray1D("b", _grid(1), _iso())
    a.stack(b, _iso(1.0))
    a.change_all_inter_connection_type(b, _iso(2.0))
    sys = ODESystem(Model(arrays=[a, b]).flatten())
    ia, ib = sys.state.index["a-0"]["A"], sys.state.index["b-0"]["A"]
    sys.set_q([ia], 100.0)
    sys.set_q([ib], 0.0)
    res = sys.run(5.0)
    assert res.y[ia][-1] + res.y[ib][-1] == pytest.approx(100.0, rel=5e-3)
    assert res.y[ia] == pytest.approx(50.0 + 50.0 * np.exp(-4.0 * res.t), rel=5e-3)


def test_to_graph_returns_compartment_graph():
    """to_graph returns a DiGraph: a node per compartment, an edge per link carrying k_out per species."""
    arr = CompartmentArray1D("m", _grid(3), _iso(0.5))
    g = arr.to_graph()
    assert isinstance(g, nx.DiGraph)
    assert set(g.nodes) == {"m-0", "m-1", "m-2"}
    assert set(g.edges) == {("m-0", "m-1"), ("m-1", "m-0"), ("m-1", "m-2"), ("m-2", "m-1")}
    assert g.edges["m-0", "m-1"]["A"] == pytest.approx(0.5)


def test_to_graph_keeps_inter_array_links():
    """Links to a stacked array appear as one-way edges to that array's compartments.
    """
    a, b, old, inter = _stacked_pair()
    g = a.to_graph()
    assert set(g.nodes) == {"a-0", "a-1", "a-2", "b-0", "b-1", "b-2"}
    assert g.number_of_edges() == 7
    assert g.edges["a-1", "b-1"]["A"] == pytest.approx(0.2)
    assert not g.has_edge("b-1", "a-1")


def test_to_graph_with_unresolved_ficks_links():
    """Unresolved FicksConnections (no rates yet) give edges without rate attributes."""
    arr = CompartmentArray1D("m", _grid(2), FicksConnection({"A": 1.0 * unit.nm ** 2 / unit.sec}))
    g = arr.to_graph()
    assert set(g.edges) == {("m-0", "m-1"), ("m-1", "m-0")}
    assert g.edges["m-0", "m-1"] == {}
