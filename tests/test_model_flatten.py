"""Tests for openrxn.model: Model array/compartment bookkeeping (periodicity,
box length, duplicate IDs), flatten() ID generation, missing-compartment
detection, resolution of FicksConnection and ResConnection into
DivByVConnection (face area and centre distance inferred from positions,
periodic wrapping), and FlatModel helpers (n_compartments, add_rxn, to_graph).
"""
import numpy as np
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment, Compartment3D, Reservoir
from openrxn.compartments.arrays import CompartmentArray1D, CompartmentArray3D
from openrxn.connections import IsotropicConnection, FicksConnection, DivByVConnection, ResConnection
from openrxn.model import Model
from openrxn.model_templates import OneCompartmentModel

NM = unit.nm
D = 1.0 * unit.nm ** 2 / unit.sec


def _box3d(ID, x, y, z):
    """Compartment3D spanning the given (lo, hi) nm intervals, with face areas set."""
    lx, ly, lz = x[1] - x[0], y[1] - y[0], z[1] - z[0]
    sa = {"yz": ly * lz * NM ** 2, "xz": lx * lz * NM ** 2, "xy": lx * ly * NM ** 2}
    pos = [(x[0] * NM, x[1] * NM), (y[0] * NM, y[1] * NM), (z[0] * NM, z[1] * NM)]
    return Compartment3D(ID, pos=pos, surface_area=sa)


def _kv(flat, c1, c2, s="A"):
    """Resolved (k_out, k_in) of the c1 -> c2 connection, in nm^3/s."""
    _, conn = flat.compartments[c1].connections[c2]
    assert isinstance(conn, DivByVConnection)
    return tuple(k.to(NM ** 3 / unit.sec).magnitude for k in conn.species_rates[s])


def _array3d(ID, x, y, z, conn, periodic=(False, False, False)):
    return CompartmentArray3D(ID, np.asarray(x) * NM, np.asarray(y) * NM, np.asarray(z) * NM,
                              conn, periodic=list(periodic))


def test_flatten_keeps_standalone_ids():
    """A standalone compartment keeps its plain ID after flattening."""
    A = Species("A")
    c = Compartment("main")
    c.add_rxn_to_compartment(Reaction("d", [A], [], [1], [], kf=1.0 / unit.sec))
    flat = Model(compartments=[c]).flatten()
    assert list(flat.compartments) == ["main"]


def test_flatten_array_ids():
    """Array compartments are flattened to '<array_ID>-<index>' IDs."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec}, dim=1)
    arr = CompartmentArray1D("m", np.linspace(0, 4, 5) * unit.nm, conn)
    flat = Model(arrays=[arr]).flatten()
    assert len(flat.compartments) == 4
    assert "m-0" in flat.compartments
    assert "m-3" in flat.compartments


def test_missing_compartment_raises():
    """flatten() raises if a connection references a compartment not in the model."""
    c1 = Compartment("c1")
    c2 = Compartment("c2")
    c1.connect(c2, IsotropicConnection({"A": 0.1 / unit.sec}))
    with pytest.raises(ValueError):
        Model(compartments=[c1]).flatten()  # c2 never added


def test_ficks_resolved_to_divbyv():
    """A FicksConnection (with area+distance) becomes a DivByVConnection on flatten."""
    sa = {"xy": 1 * unit.nm ** 2, "yz": 1 * unit.nm ** 2, "xz": 1 * unit.nm ** 2}
    c1 = Compartment3D("c1", pos=[(0 * unit.nm, 1 * unit.nm),
                                  (0 * unit.nm, 1 * unit.nm),
                                  (0 * unit.nm, 1 * unit.nm)], surface_area=sa)
    c2 = Compartment3D("c2", pos=[(1 * unit.nm, 2 * unit.nm),
                                  (0 * unit.nm, 1 * unit.nm),
                                  (0 * unit.nm, 1 * unit.nm)], surface_area=sa)
    conn = FicksConnection({"A": 2 * unit.nm ** 2 / unit.sec},
                           surface_area=3 * unit.nm ** 2, ic_distance=6 * unit.nm)
    c1.connect(c2, conn)
    c2.connect(c1, conn)
    flat = Model(compartments=[c1, c2]).flatten()
    _, resolved = flat.compartments["c1"].connections["c2"]
    assert isinstance(resolved, DivByVConnection)
    assert resolved.species_rates["A"][0].to(unit.nm ** 3 / unit.sec).magnitude == pytest.approx(1.0)


def test_one_compartment_model_initializes_and_flattens():
    """OneCompartmentModel imports, sets up base Model attrs, and flattens."""
    m = OneCompartmentModel()
    assert "main" in m.compartments
    assert m.arrays == {}        # proves super().__init__ ran
    assert m.periodic is None
    assert "main" in m.flatten().compartments


# --- Model.add_array / add_compartment ---------------------------------------

def test_add_array_stacks_box_length_along_non_periodic_dim():
    """Arrays stacked along a non-periodic z add their z lengths; periodic x, y stay fixed."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec})
    lower = _array3d("lower", [0, 1, 2], [0, 1, 2], [0, 1], conn, periodic=(True, True, False))
    upper = _array3d("upper", [0, 1, 2], [0, 1, 2], [1, 4], conn, periodic=(True, True, False))
    m = Model(arrays=[lower, upper])
    assert m.periodic == [True, True, False]
    assert [b.to(NM).magnitude for b in m.box_len] == pytest.approx([2.0, 2.0, 4.0])
    assert set(m.arrays) == {"lower", "upper"}


def test_add_array_rejects_mismatched_periodicity():
    """Adding an array whose periodicity differs from the model's raises ValueError."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec})
    a = _array3d("a", [0, 1], [0, 1], [0, 1], conn, periodic=(True, False, False))
    b = _array3d("b", [0, 1], [0, 1], [0, 1], conn, periodic=(False, False, False))
    with pytest.raises(ValueError, match="periodicity"):
        Model(arrays=[a, b])


def test_add_array_rejects_different_length_along_periodic_dim():
    """Arrays of different x length cannot share a periodic x dimension."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec})
    a = _array3d("a", [0, 1], [0, 1], [0, 1], conn, periodic=(True, False, False))
    b = _array3d("b", [0, 2], [0, 1], [1, 2], conn, periodic=(True, False, False))
    with pytest.raises(ValueError, match="periodic dimension"):
        Model(arrays=[a, b])


def test_add_array_rejects_duplicate_id():
    """Two arrays with the same array_ID raise ValueError."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec})
    a = _array3d("same", [0, 1], [0, 1], [0, 1], conn)
    b = _array3d("same", [0, 1], [0, 1], [0, 1], conn)
    with pytest.raises(ValueError, match="Duplicate array ID"):
        Model(arrays=[a, b])


def test_add_compartment_rejects_duplicate_id():
    """Two standalone compartments with the same ID raise ValueError."""
    with pytest.raises(ValueError, match="Duplicate compartment ID"):
        Model(compartments=[Compartment("c"), Compartment("c")])


def test_add_array_stacks_two_1d_arrays():
    """Two non-periodic 1D arrays of different length stack: box length 4 + 2 nm."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec}, dim=1)
    a = CompartmentArray1D("a", np.linspace(0, 4, 5) * NM, conn)
    b = CompartmentArray1D("b", np.linspace(4, 6, 3) * NM, conn)
    m = Model(arrays=[a, b])
    assert m.box_len[0].to(NM).magnitude == pytest.approx(6.0)


def test_add_array_does_not_mutate_input_array():
    """Building a model leaves the arrays' own box_len unchanged, so they can be reused."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec})
    lower = _array3d("lower", [0, 1], [0, 1], [0, 1], conn, periodic=(True, True, False))
    upper = _array3d("upper", [0, 1], [0, 1], [1, 3], conn, periodic=(True, True, False))
    Model(arrays=[lower, upper])
    assert lower.box_len[2].to(NM).magnitude == pytest.approx(1.0)
    assert Model(arrays=[lower]).box_len[2].to(NM).magnitude == pytest.approx(1.0)


# --- Model.resolve_ficks -------------------------------------------------------

@pytest.mark.parametrize("axis, c2_span, area", [
    ("x", ([1, 2], [0, 4], [0, 6]), 6.0),
    ("y", ([0, 2], [2, 4], [0, 6]), 3.0),
    ("z", ([0, 2], [0, 4], [3, 6]), 2.0),
])
def test_ficks_area_from_adjoining_face(axis, c2_span, area):
    """With no surface_area, the smaller of the two shared faces (yz, xz or xy) is used."""
    c1 = _box3d("c1", [0, 1], [0, 2], [0, 3])     # faces: yz=6, xz=3, xy=2 nm^2
    c2 = _box3d("c2", *c2_span)                    # larger faces, so min() picks c1's
    conn = FicksConnection({"A": D}, ic_distance=1.0 * NM)
    c1.connect(c2, conn)
    c2.connect(c1, conn)
    flat = Model(compartments=[c1, c2]).flatten()
    assert _kv(flat, "c1", "c2") == pytest.approx((area, area))
    assert _kv(flat, "c2", "c1") == pytest.approx((area, area))


def test_ficks_non_adjoining_compartments_raise():
    """Ficks resolution without surface_area fails for compartments that share no face."""
    c1 = _box3d("c1", [0, 1], [0, 1], [0, 1])
    c2 = _box3d("c2", [5, 6], [5, 6], [5, 6])
    conn = FicksConnection({"A": D}, ic_distance=1.0 * NM)
    c1.connect(c2, conn)
    c2.connect(c1, conn)
    with pytest.raises(ValueError, match="adjoining face"):
        Model(compartments=[c1, c2]).flatten()


def test_ficks_overlapping_compartments_raise():
    """Two different compartments at the same position have no centre distance: ValueError."""
    c1 = _box3d("c1", [0, 1], [0, 1], [0, 1])
    c2 = _box3d("c2", [0, 1], [0, 1], [0, 1])
    conn = FicksConnection({"A": D}, surface_area=1.0 * NM ** 2)
    c1.connect(c2, conn)
    c2.connect(c1, conn)
    with pytest.raises(ValueError, match="same position"):
        Model(compartments=[c1, c2]).flatten()


def test_ficks_area_and_distance_from_array_positions():
    """In a 2x1x1 array of 1x2x3 nm cells, Fick's rate is D * (yz = 6 nm^2) / (dx = 1 nm)."""
    arr = _array3d("g", [0, 1, 2], [0, 2], [0, 3], FicksConnection({"A": D}))
    flat = Model(arrays=[arr]).flatten()
    assert _kv(flat, "g-0_0_0", "g-1_0_0") == pytest.approx((6.0, 6.0))
    assert _kv(flat, "g-1_0_0", "g-0_0_0") == pytest.approx((6.0, 6.0))


def test_ficks_distance_wraps_across_periodic_boundary():
    """Cells 0 and 2 of a periodic 3-cell x ring are 1 nm apart (not 2), in both directions."""
    arr = _array3d("p", [0, 1, 2, 3], [0, 2], [0, 3], IsotropicConnection({"A": 0.1 / unit.sec}),
                   periodic=(True, False, False))
    c0, c2 = arr.compartments[(0, 0, 0)], arr.compartments[(2, 0, 0)]
    fwd = FicksConnection({"A": D}, surface_area=6.0 * NM ** 2)
    bwd = FicksConnection({"A": D}, surface_area=6.0 * NM ** 2)
    c0.connect(c2, fwd, warn_overwrite=False)
    c2.connect(c0, bwd, warn_overwrite=False)
    flat = Model(arrays=[arr]).flatten()
    assert _kv(flat, "p-0_0_0", "p-2_0_0") == pytest.approx((6.0, 6.0))
    assert _kv(flat, "p-2_0_0", "p-0_0_0") == pytest.approx((6.0, 6.0))


def test_ficks_periodic_ring_resolves_every_pair():
    """Every link of a periodic 3-cell x ring of 1x2x3 nm cells, wrap included, is D * 6 nm^2 / 1 nm."""
    arr = _array3d("p", [0, 1, 2, 3], [0, 2], [0, 3], FicksConnection({"A": D}), periodic=(True, False, False))
    flat = Model(arrays=[arr]).flatten()
    links = [(c, n) for c in flat.compartments for n in flat.compartments[c].connections]
    assert len(links) == 6
    for c, n in links:
        assert _kv(flat, c, n) == pytest.approx((6.0, 6.0))


def test_ficks_single_cell_periodic_axis_links_to_own_image():
    """A one-cell periodic x axis links each cell to its image one box length away (D * yz / Lx)."""
    arr = _array3d("col", [0, 2], [0, 2], [0, 1, 2], FicksConnection({"A": D}), periodic=(True, True, False))
    flat = Model(arrays=[arr]).flatten()
    assert _kv(flat, "col-0_0_0", "col-0_0_1") == pytest.approx((4.0, 4.0))   # xy = 4 nm^2, dz = 1 nm
    assert _kv(flat, "col-0_0_0", "col-0_0_0") == pytest.approx((1.0, 1.0))   # yz = 2 nm^2, Lx = 2 nm


def test_ficks_rates_per_face_in_non_cubic_array():
    """In a 2x1x2 array of 1x1x2 nm cells, x-neighbours get D*2/1 and z-neighbours D*1/2."""
    arr = _array3d("g", [0, 1, 2], [0, 1], [0, 2, 4], FicksConnection({"A": D}))
    flat = Model(arrays=[arr]).flatten()
    assert _kv(flat, "g-0_0_0", "g-1_0_0") == pytest.approx((2.0, 2.0))
    assert _kv(flat, "g-0_0_0", "g-0_0_1") == pytest.approx((0.5, 0.5))


def test_ficks_area_across_periodic_boundary():
    """A fresh Ficks connection between the end cells of a periodic x ring uses the yz face."""
    arr = _array3d("p", [0, 1, 2, 3], [0, 2], [0, 3], IsotropicConnection({"A": 0.1 / unit.sec}),
                   periodic=(True, False, False))
    c0, c2 = arr.compartments[(0, 0, 0)], arr.compartments[(2, 0, 0)]
    conn = FicksConnection({"A": D}, ic_distance=1.0 * NM)
    c0.connect(c2, conn, warn_overwrite=False)
    c2.connect(c0, conn, warn_overwrite=False)
    flat = Model(arrays=[arr]).flatten()
    assert _kv(flat, "p-0_0_0", "p-2_0_0") == pytest.approx((6.0, 6.0))


def test_ficks_distance_between_standalone_compartments():
    """Two standalone 3D boxes (centres 0.5 and 2.0 nm) resolve to D * 1 nm^2 / 1.5 nm."""
    c1 = _box3d("c1", [0, 1], [0, 1], [0, 1])
    c2 = _box3d("c2", [1, 3], [0, 1], [0, 1])
    conn = FicksConnection({"A": D})
    c1.connect(c2, conn)
    c2.connect(c1, conn)
    flat = Model(compartments=[c1, c2]).flatten()
    assert _kv(flat, "c1", "c2") == pytest.approx((1 / 1.5, 1 / 1.5))


# --- Model.resolve_res -----------------------------------------------------------

@pytest.mark.parametrize("face, kv", [("x", 6.0 / 0.5), ("y", 3.0 / 1.0), ("z", 2.0 / 1.5)])
def test_res_connection_geometry_from_face(face, kv):
    """ResConnection(face) uses that face's area and half the cell width along it: kV = D*A/(width/2)."""
    c = _box3d("c", [0, 1], [0, 2], [0, 3])         # yz=6, xz=3, xy=2 nm^2
    res = Reservoir("res", conc_funcs={"A": lambda t: 0.0})
    c.connect(res, ResConnection({"A": D}, face=face))
    flat = Model(compartments=[c, res]).flatten()
    neighbour, _ = flat.compartments["c"].connections["res"]
    assert isinstance(neighbour, Reservoir)
    assert _kv(flat, "c", "res") == pytest.approx((kv, kv))


def test_shared_res_connection_resolved_per_compartment():
    """One ResConnection(face='z') shared by cells 1 and 2 nm tall gives each its own D*A/(width/2)."""
    lo = _box3d("lo", [0, 1], [0, 1], [0, 1])
    hi = _box3d("hi", [1, 2], [0, 1], [0, 2])
    res = Reservoir("res", conc_funcs={"A": lambda t: 0.0})
    conn = ResConnection({"A": D}, face="z")
    lo.connect(res, conn)
    hi.connect(res, conn)
    flat = Model(compartments=[lo, hi, res]).flatten()
    assert _kv(flat, "lo", "res") == pytest.approx((2.0, 2.0))
    assert _kv(flat, "hi", "res") == pytest.approx((1.0, 1.0))


def test_res_connection_without_face_or_area_raises():
    """Without face or surface_area the reservoir interface cannot be resolved."""
    c = _box3d("c", [0, 1], [0, 1], [0, 1])
    res = Reservoir("res", conc_funcs={"A": lambda t: 0.0})
    c.connect(res, ResConnection({"A": D}, ic_distance=1.0 * NM))
    with pytest.raises(ValueError, match="surface area"):
        Model(compartments=[c, res]).flatten()


def test_res_connection_without_face_or_distance_raises():
    """With surface_area but no face or ic_distance the distance cannot be resolved."""
    c = _box3d("c", [0, 1], [0, 1], [0, 1])
    res = Reservoir("res", conc_funcs={"A": lambda t: 0.0})
    c.connect(res, ResConnection({"A": D}, surface_area=1.0 * NM ** 2))
    with pytest.raises(ValueError, match="intercompartmental distance"):
        Model(compartments=[c, res]).flatten()


# --- FlatModel -----------------------------------------------------------------------

def test_flat_id_collision_raises():
    """A standalone 'm-0' collides with compartment 0 of array 'm' after flattening."""
    arr = CompartmentArray1D("m", np.linspace(0, 2, 3) * NM, IsotropicConnection({"A": 0.1 / unit.sec}, dim=1))
    with pytest.raises(ValueError, match="m-0"):
        Model(arrays=[arr], compartments=[Compartment("m-0")]).flatten()


def test_flat_n_compartments_counts_arrays_and_standalone():
    """n_compartments counts array cells plus standalone compartments."""
    arr = CompartmentArray1D("m", np.linspace(0, 4, 5) * NM, IsotropicConnection({"A": 0.1 / unit.sec}, dim=1))
    flat = Model(arrays=[arr], compartments=[Compartment("extra")]).flatten()
    assert flat.n_compartments() == 5


def test_flat_add_rxn_to_selected_and_all():
    """add_rxn adds to the listed compartments only, or to every one with 'all'."""
    A, B = Species("A"), Species("B")
    r1 = Reaction("r1", [A], [B], [1], [1], kf=1.0 / unit.sec)
    r2 = Reaction("r2", [B], [A], [1], [1], kf=1.0 / unit.sec)
    arr = CompartmentArray1D("m", np.linspace(0, 3, 4) * NM, IsotropicConnection({"A": 0.1 / unit.sec}, dim=1))
    flat = Model(arrays=[arr]).flatten()
    flat.add_rxn(r1, ["m-1"])
    assert [[r.ID for r in flat.compartments[c].reactions] for c in ["m-0", "m-1", "m-2"]] == [[], ["r1"], []]
    flat.add_rxn(r2)
    assert [[r.ID for r in flat.compartments[c].reactions] for c in ["m-0", "m-1", "m-2"]] == \
        [["r2"], ["r1", "r2"], ["r2"]]
    with pytest.raises(AssertionError, match="not in Model"):
        flat.add_rxn(r1, ["nope"])


def test_to_graph_nodes_positions_and_edges():
    """to_graph: one node per compartment at the projected centre, one edge per connection with k_out."""
    arr = _array3d("g", [0, 2, 4], [0, 4], [0, 6], IsotropicConnection({"A": 0.4 / unit.sec}))
    flat = Model(arrays=[arr]).flatten()
    g = flat.to_graph(scale=10)
    assert set(g.nodes) == {"g-0_0_0", "g-1_0_0"}
    # centres (1,2,3) and (3,2,3) nm, x10; projection (x - 0.7 y, z + 1.2 y)
    assert g.nodes["g-0_0_0"]["viz"]["position"] == pytest.approx({"x": -4.0, "y": 54.0})
    assert g.nodes["g-1_0_0"]["viz"]["position"] == pytest.approx({"x": 16.0, "y": 54.0})
    assert set(g.edges) == {("g-0_0_0", "g-1_0_0"), ("g-1_0_0", "g-0_0_0")}
    assert g.edges["g-0_0_0", "g-1_0_0"]["A"] == pytest.approx(0.4)


def test_to_graph_1d_array():
    """to_graph of a 1D array has one node per cell (y = z = 0) and edges between neighbours."""
    arr = CompartmentArray1D("m", np.linspace(0, 3, 4) * NM, IsotropicConnection({"A": 0.1 / unit.sec}, dim=1))
    g = Model(arrays=[arr]).flatten().to_graph()
    assert set(g.nodes) == {"m-0", "m-1", "m-2"}
    assert set(g.edges) == {("m-0", "m-1"), ("m-1", "m-0"), ("m-1", "m-2"), ("m-2", "m-1")}
    assert g.nodes["m-2"]["viz"]["position"] == pytest.approx({"x": 25.0, "y": 0.0})


def test_to_graph_places_compartments_without_position_at_origin():
    """Reservoirs and plain compartments have no pos: their nodes sit at (0, 0) and their edges are kept."""
    c = _box3d("c", [0, 2], [0, 1], [0, 1])
    res = Reservoir("res", conc_funcs={"A": lambda t: 0.0})
    plain = Compartment("plain", volume=1.0 * NM ** 3)
    c.connect(res, ResConnection({"A": D}, face="x"))
    plain.connect(c, IsotropicConnection({"A": 0.3 / unit.sec}))
    c.connect(plain, IsotropicConnection({"A": 0.3 / unit.sec}))
    g = Model(compartments=[c, res, plain]).flatten().to_graph()
    for node in ["res", "plain"]:
        assert g.nodes[node]["viz"]["position"] == {"x": 0.0, "y": 0.0}
    assert g.edges["c", "res"]["A"] == pytest.approx(1.0)      # D * 1 nm^2 / (2 nm / 2)
    assert g.edges["plain", "c"]["A"] == pytest.approx(0.3)
