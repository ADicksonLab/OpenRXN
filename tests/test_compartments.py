"""Tests for openrxn.compartments: volume computation and position validation,
reaction handling, connection add/overwrite/remove, copy semantics, repr,
array topology, and Reservoir concentration sources and copies.

"""
import logging

import numpy as np
import pytest
from pint import DimensionalityError

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment, Compartment1D, Compartment2D, Compartment3D, Reservoir
from openrxn.compartments.arrays import CompartmentArray1D, CompartmentArray2D, CompartmentArray3D
from openrxn.connections import IsotropicConnection


def test_compartment_1d_volume():
    """Compartment1D volume is the length of its single position interval."""
    c = Compartment1D("c", pos=[(0 * unit.nm, 2 * unit.nm)])
    assert c.volume.to(unit.nm).magnitude == pytest.approx(2.0)


def test_compartment_2d_volume():
    """Compartment2D volume is the product of its two interval lengths."""
    c = Compartment2D("c", pos=[(0 * unit.nm, 2 * unit.nm), (0 * unit.nm, 3 * unit.nm)])
    assert c.volume.to(unit.nm ** 2).magnitude == pytest.approx(6.0)


def test_compartment_3d_volume():
    """Compartment3D volume is the product of its three interval lengths."""
    c = Compartment3D("c", pos=[(0 * unit.nm, 2 * unit.nm),
                                 (0 * unit.nm, 3 * unit.nm),
                                 (0 * unit.nm, 4 * unit.nm)])
    assert c.volume.to(unit.nm ** 3).magnitude == pytest.approx(24.0)


def test_duplicate_reaction_ignored():
    """Adding the same reaction twice leaves the compartment with a single copy."""
    A = Species("A")
    r = Reaction("decay", [A], [], [1], [], kf=1.0 / unit.sec)
    c = Compartment("c")
    c.add_rxn_to_compartment(r)
    c.add_rxn_to_compartment(r)
    assert len(c.reactions) == 1


def test_copy_preserves_reactions_new_id():
    """copy(ID=...) returns a compartment with a new ID but the same reactions."""
    A = Species("A")
    r = Reaction("decay", [A], [], [1], [], kf=1.0 / unit.sec)
    c = Compartment("c")
    c.add_rxn_to_compartment(r)
    c2 = c.copy(ID="c2")
    assert c2.ID == "c2"
    assert len(c2.reactions) == 1


def test_array_1d_size_and_topology():
    """A 1D array of N intervals has N compartments; ends have 1 neighbor, interior 2."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec}, dim=1)
    arr = CompartmentArray1D("m", np.linspace(0, 4, 5) * unit.nm, conn)
    assert len(arr.compartments) == 4
    assert len(arr.compartments[0].connections) == 1
    assert len(arr.compartments[1].connections) == 2


def test_array_1d_periodic_wraps():
    """With periodic=True the end compartments gain the wrap-around neighbor."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec}, dim=1)
    arr = CompartmentArray1D("m", np.linspace(0, 4, 5) * unit.nm, conn, periodic=True)
    assert len(arr.compartments[0].connections) == 2
    assert len(arr.compartments[3].connections) == 2


def test_array_3d_size():
    """A 2x2x2 grid has 8 compartments."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec})
    arr = CompartmentArray3D(
        "g",
        np.linspace(0, 2, 3) * unit.nm,
        np.linspace(0, 2, 3) * unit.nm,
        np.linspace(0, 2, 3) * unit.nm,
        conn,
    )
    assert len(arr.compartments) == 8


def test_array_2d_grid_size():
    """A 3x3 mesh should contain nx*ny = 9 compartments."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec}, dim=2)
    arr = CompartmentArray2D(
        "g", np.linspace(0, 3, 4) * unit.nm, np.linspace(0, 3, 4) * unit.nm, conn
    )
    assert len(arr.compartments) == arr.nx * arr.ny == 9


def test_reservoir_time_conc():
    """A time-dependent conc_funcs entry is stored and callable."""
    res = Reservoir("r", conc_funcs={"A": lambda t: 1.0 + t})
    assert res.conc_funcs["A"](0) == pytest.approx(1.0)
    assert res.conc_funcs["A"](2) == pytest.approx(3.0)


def test_reservoir_static_conc_per_key():
    """Each static concentration should be retrievable through its own key."""
    res = Reservoir("r", concs={"A": 2.0 / unit.nm ** 3, "B": 5.0 / unit.nm ** 3})
    assert res.conc_funcs["A"](0) == pytest.approx(2.0)
    assert res.conc_funcs["B"](0) == pytest.approx(5.0)

def _conn(k):
    return IsotropicConnection({"A": k / unit.sec})


def _two_rxns():
    A, B = Species("A"), Species("B")
    r1 = Reaction("r1", [A], [B], [1], [1], kf=1.0 / unit.sec)
    r2 = Reaction("r2", [B], [], [1], [], kf=2.0 / unit.sec)
    return r1, r2


def test_add_rxns_keeps_order_and_warns_on_repeat(caplog):
    """add_rxns_to_compartment keeps insertion order and warns when a reaction ID repeats."""
    r1, r2 = _two_rxns()
    c = Compartment("c")
    with caplog.at_level(logging.WARNING):
        c.add_rxns_to_compartment([r1, r2, r1])
    assert [r.ID for r in c.reactions] == ["r1", "r2"]
    assert "Reaction r1 already in compartment c" in caplog.text


def test_show_all_rxns():
    """show_all_rxns returns one display string per reaction, in insertion order."""
    r1, r2 = _two_rxns()
    c = Compartment("c")
    assert c.show_all_rxns() == []
    c.add_rxns_to_compartment([r1, r2])
    assert c.show_all_rxns() == [r1.display(), r2.display()]
    assert c.show_all_rxns()[0].startswith("A ---> B")


def test_connect_overwrite_warns_and_replaces(caplog):
    """Re-connecting a neighbor replaces the link and logs a warning naming both array tags."""
    c = Compartment(0, array_ID="m")
    n = Compartment(1, array_ID="m")
    old, new = _conn(0.1), _conn(0.2)
    c.connect(n, old)
    with caplog.at_level(logging.WARNING):
        c.connect(n, new)
    assert c.connections == {"m-1": (n, new)}
    assert "Overwriting connection between m-0 and m-1" in caplog.text


def test_connect_overwrite_silent_when_disabled(caplog):
    """warn_overwrite=False replaces the link without logging a warning."""
    c, n = Compartment("c"), Compartment("n")
    new = _conn(0.2)
    c.connect(n, _conn(0.1))
    with caplog.at_level(logging.WARNING):
        c.connect(n, new, warn_overwrite=False)
    assert c.connections == {"n": (n, new)}
    assert "Overwriting" not in caplog.text


def test_remove_connection_drops_only_that_neighbor():
    """remove_connection deletes the link to one neighbor and keeps the others."""
    c, n1, n2 = Compartment("c"), Compartment("n1"), Compartment("n2")
    conn = _conn(0.1)
    c.connect(n1, conn)
    c.connect(n2, conn)
    c.remove_connection(n1)
    assert c.connections == {"n2": (n2, conn)}


def test_remove_missing_connection_only_warns(caplog):
    """Removing a link that does not exist logs a warning and leaves the links unchanged."""
    c, n1 = Compartment("c"), Compartment("n1")
    conn = _conn(0.1)
    c.connect(n1, conn)
    with caplog.at_level(logging.WARNING):
        c.remove_connection(Compartment("ghost"))
    assert "does not exist" in caplog.text
    assert c.connections == {"n1": (n1, conn)}


def test_copy_defaults_and_independence():
    """copy() keeps ID, array_ID, volume and links; the copy's containers are independent."""
    r1, r2 = _two_rxns()
    c = Compartment("c", array_ID="arr", volume=2.0 * unit.L)
    n = Compartment("n")
    c.add_rxn_to_compartment(r1)
    c.connect(n, _conn(0.1))
    cp = c.copy()
    assert (cp.ID, cp.array_ID) == ("c", "arr")
    assert cp.volume.to(unit.L).magnitude == pytest.approx(2.0)
    assert cp.connections == c.connections
    cp.add_rxn_to_compartment(r1)
    assert [r.ID for r in cp.reactions] == ["r1"]
    cp.add_rxn_to_compartment(r2)
    cp.remove_connection(n)
    assert [r.ID for r in c.reactions] == ["r1"]
    assert "n" in c.connections
    assert c.copy(delete_array_ID=True).array_ID is None


def test_copy_3d_keeps_surface_area():
    """Copying a Compartment3D keeps its class, surface areas, position and volume."""
    sa = {"xy": 6 * unit.nm ** 2, "yz": 12 * unit.nm ** 2, "xz": 8 * unit.nm ** 2}
    c = Compartment3D("c", pos=[(0 * unit.nm, 2 * unit.nm),
                                (0 * unit.nm, 3 * unit.nm),
                                (0 * unit.nm, 4 * unit.nm)], surface_area=sa)
    cp = c.copy(ID="c2")
    assert isinstance(cp, Compartment3D)
    assert cp.ID == "c2"
    assert cp.surface_area == sa
    assert cp.pos == c.pos
    assert cp.volume.to(unit.nm ** 3).magnitude == pytest.approx(24.0)


def test_repr_bare_compartment():
    """repr of an empty compartment shows ID, volume, no position and zero counts."""
    c = Compartment("c", volume=2 * unit.L)
    assert repr(c) == "Compartment(ID='c', volume=2 liter, pos=None, n_rxns=0, n_connections=0)"


def test_repr_counts_reactions_and_connections():
    """repr reports the position bounds and the numbers of reactions and connections."""
    r1, r2 = _two_rxns()
    c = Compartment1D(0, pos=[(0 * unit.nm, 2 * unit.nm)])
    c.add_rxns_to_compartment([r1, r2])
    c.connect(Compartment("n"), _conn(0.1))
    s = repr(c)
    assert s.startswith("Compartment(ID=0, volume=2 nanometer, pos=((")
    assert s.endswith("n_rxns=2, n_connections=1)")


@pytest.mark.parametrize("cls, n_dims", [(Compartment1D, 2), (Compartment2D, 1), (Compartment3D, 2)])
def test_spatial_compartment_rejects_wrong_dimension(cls, n_dims):
    """Compartment1D/2D/3D require exactly 1/2/3 position intervals."""
    with pytest.raises(AssertionError):
        cls("c", pos=[(0 * unit.nm, 1 * unit.nm)] * n_dims)


@pytest.mark.parametrize("cls, n_dims", [(Compartment1D, 1), (Compartment2D, 2), (Compartment3D, 3)])
def test_spatial_compartment_requires_length_units(cls, n_dims):
    """Positions in non-length units are rejected when the volume is checked against nm**d."""
    with pytest.raises(DimensionalityError):
        cls("c", pos=[(0 * unit.sec, 1 * unit.sec)] * n_dims)


def test_reservoir_static_conc_converted_to_per_nm3():
    """Static concentrations become plain floats in molecules per nm**3, constant in time."""
    res = Reservoir("r", concs={"A": 1e9 / unit.um ** 3})
    val = res.conc_funcs["A"](0.0)
    assert isinstance(val, float)
    assert val == pytest.approx(1.0)
    assert res.conc_funcs["A"](100.0) == val


def test_reservoir_user_conc_func_is_kept_unchanged():
    """A user-supplied conc_funcs entry is stored as the same callable, not wrapped."""
    f = lambda t: 50.0 * unit.nm ** -3
    res = Reservoir("r", conc_funcs={"A": f})
    assert res.conc_funcs["A"] is f


def test_reservoir_rejects_key_in_both_sources():
    """A species given in both concs and conc_funcs raises ValueError."""
    with pytest.raises(ValueError, match="both concs and conc_funcs"):
        Reservoir("r", concs={"A": 1.0 / unit.nm ** 3}, conc_funcs={"A": lambda t: 1.0})


def test_reservoir_rejects_non_callable_conc_func():
    """A conc_funcs entry that is not callable raises ValueError."""
    with pytest.raises(ValueError, match="not callable"):
        Reservoir("r", conc_funcs={"A": 1.0 / unit.nm ** 3})


def test_reservoir_copy_keeps_concs_and_ids():
    """Reservoir.copy keeps ID, array_ID and concentration functions in a new dict."""
    f = lambda t: 2.0 * t
    res = Reservoir("r", array_ID="bath", concs={"A": 3.0 / unit.nm ** 3}, conc_funcs={"B": f})
    cp = res.copy()
    assert isinstance(cp, Reservoir)
    assert (cp.ID, cp.array_ID) == ("r", "bath")
    assert cp.conc_funcs["A"](5.0) == pytest.approx(3.0)
    assert cp.conc_funcs["B"] is f
    cp2 = res.copy(ID="r2", delete_array_ID=True)
    assert (cp2.ID, cp2.array_ID) == ("r2", None)
    cp2.conc_funcs["C"] = f
    assert set(res.conc_funcs) == {"A", "B"}


def test_reservoir_copy_keeps_pos_and_volume():
    """Reservoir.copy, like Compartment.copy, preserves position and volume."""
    res = Reservoir("r", pos=[(0 * unit.nm, 5 * unit.nm)], volume=3.0 * unit.L,
                    concs={"A": 1.0 / unit.nm ** 3})
    cp = res.copy()
    assert cp.pos == res.pos
    assert cp.volume is not None
    assert cp.volume.to(unit.L).magnitude == pytest.approx(3.0)


def test_reservoir_copy_drops_links_and_reactions():
    """Reservoir.copy leaves out connections and reactions, so a flattened reservoir stays a fixed source."""
    r1, _ = _two_rxns()
    res = Reservoir("r", concs={"A": 1.0 / unit.nm ** 3})
    res.add_rxn_to_compartment(r1)
    res.connect(Compartment("c"), _conn(0.1))
    cp = res.copy()
    assert cp.connections == {}
    assert cp.reactions == []
