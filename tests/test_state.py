"""Tests for openrxn.systems.state.State: the species index, sparse
inclusion of only participating species, and compartment position arrays.

Also covers building a State from a dataframe and exporting it with
to_dataframe / to_csv / to_csv_no_q (files go to pytest's tmp_path).

NOTE: the per-compartment index is built by iterating a Python set, so the
*order* of state.species / q_val is not contractual -- tests address species
through state.index[comp][species], never by assumed position.
"""
import numpy as np
import pandas as pd
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction
from openrxn.compartments.compartment import Compartment
from openrxn.compartments.arrays import CompartmentArray1D
from openrxn.connections import IsotropicConnection
from openrxn.model import Model
from openrxn.systems.state import State


def test_index_maps_species(decay_model):
    """A single-species decay model yields one state entry indexed by [comp][species]."""
    st = State(model=decay_model["flat"])
    assert "main" in st.index
    assert "A" in st.index["main"]
    assert st.size == 1
    assert st.q_val.shape == (1,)


def test_index_includes_all_reacting_species(equilibrium_model):
    """Both species of A<=>B are indexed and the state has size 2."""
    st = State(model=equilibrium_model["flat"])
    assert set(st.index["main"].keys()) == {"A", "B"}
    assert st.size == 2


def test_unused_species_excluded():
    """A species that participates in no reaction or connection is not in the state."""
    A = Species("A")
    Species("B")  # defined but never used in a reaction/connection
    c = Compartment("main")
    c.add_rxn_to_compartment(Reaction("d", [A], [], [1], [], kf=1.0 / unit.sec))
    flat = Model(compartments=[c]).flatten()
    st = State(model=flat)
    assert "A" in st.index["main"]
    assert "B" not in st.index["main"]
    assert list(st.species) == ["A"]


def test_position_arrays_from_array_model():
    """Compartment-center x positions of a 1D mesh are the interval midpoints."""
    conn = IsotropicConnection({"A": 0.1 / unit.sec}, dim=1)
    arr = CompartmentArray1D("m", np.linspace(0, 4, 5) * unit.nm, conn)
    flat = Model(arrays=[arr]).flatten()
    st = State(model=flat)
    assert st.size == 4
    assert sorted(st.x_pos) == pytest.approx([0.5, 1.5, 2.5, 3.5])

@pytest.fixture
def box_model():
    """One 3D-positioned compartment with A --> B; centre at (1, 2, 13) nm."""
    nm = unit.nm
    A, B = Species("A"), Species("B")
    c = Compartment("box", pos=[[0 * nm, 2 * nm], [0 * nm, 4 * nm], [10 * nm, 16 * nm]])
    c.add_rxn_to_compartment(Reaction("ab", [A], [B], [1], [1], kf=1.0 / unit.sec))
    return Model(compartments=[c]).flatten()


@pytest.fixture
def df_two_comps():
    """A minimal dataframe describing A in two compartments, with x positions only."""
    return pd.DataFrame({"species": ["A", "B", "A"],
                         "compartment": ["left", "left", "right"],
                         "x_pos": [0.5, 0.5, 1.5]})


def test_position_arrays_xyz(box_model):
    """x, y and z positions are the per-axis midpoints of the compartment bounds (in nm)."""
    st = State(model=box_model)
    assert st.x_pos == pytest.approx([1.0, 1.0])
    assert st.y_pos == pytest.approx([2.0, 2.0])
    assert st.z_pos == pytest.approx([13.0, 13.0])


def test_position_units_are_respected(box_model):
    """Positions are expressed in the units passed to State (angstrom = 10x nm)."""
    st = State(model=box_model, units=[unit.angstrom] * 3)
    assert st.x_pos == pytest.approx([10.0, 10.0])
    assert st.z_pos == pytest.approx([130.0, 130.0])


def test_init_from_dataframe(df_two_comps):
    """A dataframe State indexes rows by [compartment][species], zeroes q_val, and keeps given positions."""
    st = State(dataframe=df_two_comps)
    assert st.index == {"left": {"A": 0, "B": 1}, "right": {"A": 2}}
    assert st.size == 3
    assert np.array_equal(st.q_val, np.zeros(3))
    assert st.x_pos == pytest.approx([0.5, 0.5, 1.5])
    assert not hasattr(st, "y_pos")
    assert not hasattr(st, "z_pos")


@pytest.mark.parametrize("missing", ["species", "compartment"])
def test_init_from_dataframe_requires_columns(df_two_comps, missing):
    """A dataframe without 'species' or 'compartment' columns raises ValueError."""
    with pytest.raises(ValueError, match="species"):
        State(dataframe=df_two_comps.drop(columns=missing))


def test_to_dataframe(box_model):
    """to_dataframe has one row per state entry with species, compartment, q_val and positions."""
    st = State(model=box_model)
    st.q_val[st.index["box"]["A"]] = 7.0
    st.q_val[st.index["box"]["B"]] = 3.0
    df = st.to_dataframe()
    assert list(df.columns) == ["species", "compartment", "q_val", "x_pos", "y_pos", "z_pos"]
    q = dict(zip(df["species"], df["q_val"]))
    assert q == {"A": 7.0, "B": 3.0}
    assert set(df["compartment"]) == {"box"}
    assert df["z_pos"].tolist() == pytest.approx([13.0, 13.0])


def test_to_dataframe_omits_absent_positions(df_two_comps):
    """A State without y/z positions exports only the columns it has."""
    df = State(dataframe=df_two_comps).to_dataframe()
    assert list(df.columns) == ["species", "compartment", "q_val", "x_pos"]


def test_to_csv(box_model, tmp_path):
    """to_csv writes q_val alongside species, compartment and positions."""
    st = State(model=box_model)
    st.q_val[st.index["box"]["B"]] = 4.5
    path = tmp_path / "state.csv"
    st.to_csv(path)
    df = pd.read_csv(path, index_col=0)
    assert list(df.columns) == ["species", "compartment", "q_val", "x_pos", "y_pos", "z_pos"]
    assert dict(zip(df["species"], df["q_val"])) == {"A": 0.0, "B": 4.5}


def test_to_csv_no_q_round_trip(two_box_model, tmp_path):
    """to_csv_no_q drops q_val, and reading the file back rebuilds the same index."""
    st = State(model=two_box_model["flat"])
    path = tmp_path / "layout.csv"
    st.to_csv_no_q(path)
    df = pd.read_csv(path, index_col=0)
    assert "q_val" not in df.columns
    assert {"species", "compartment", "x_pos"} <= set(df.columns)
    st2 = State(dataframe=df)
    assert st2.index == st.index
    assert list(st2.species) == list(st.species)
    assert list(st2.compartment) == list(st.compartment)


def test_init_from_filtered_dataframe(df_two_comps):
    """A row subset of a dataframe (non-zero-based index) still gives a positional index."""
    sub = df_two_comps[df_two_comps["compartment"] == "right"]
    st = State(dataframe=sub)
    assert st.index == {"right": {"A": 0}}
    assert st.size == 1


def test_init_from_reordered_dataframe(df_two_comps):
    """After sorting the dataframe, index[comp][species] still points at that row of the state arrays."""
    st = State(dataframe=df_two_comps.sort_values("compartment", ascending=False))
    for i in range(st.size):
        assert st.index[st.compartment[i]][st.species[i]] == i
