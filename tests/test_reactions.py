"""Tests for openrxn.reactions: Species and Reaction construction,
validation, and the unit-helper methods.
"""
import logging
import pytest

from openrxn import unit
from openrxn.reactions import Species, Reaction


def test_species_stores_id_and_kwargs():
    """Species keeps its ID and attaches any extra kwargs as attributes."""
    s = Species("A", charge=-1, label="anion")
    assert s.ID == "A"
    assert s.charge == -1
    assert s.label == "anion"


def test_reaction_stores_fields():
    """Reaction stores its reactants/products/stoich/rates and caches the species ID lists."""
    A, B, C = Species("A"), Species("B"), Species("C")
    r = Reaction("r", [A, B], [C], [1, 2], [1], kf=3.0 / unit.sec)
    assert r.reactants == [A, B]
    assert r.products == [C]
    assert r.stoich_r == [1, 2]
    assert r.stoich_p == [1]
    assert r.reactant_IDs == ["A", "B"]
    assert r.product_IDs == ["C"]


def test_reactant_stoich_length_mismatch():
    """A reactant list that doesn't match its stoichiometry list raises ValueError."""
    A, B = Species("A"), Species("B")
    with pytest.raises(ValueError):
        Reaction("r", [A, B], [], [1], [], kf=1.0 / unit.sec)


def test_product_stoich_length_mismatch():
    """A product list that doesn't match its stoichiometry list raises ValueError."""
    A, B = Species("A"), Species("B")
    with pytest.raises(ValueError):
        Reaction("r", [A], [B], [1], [1, 1], kf=1.0 / unit.sec)


def test_non_species_reactant_raises():
    """Passing a non-Species object as a reactant raises TypeError."""
    with pytest.raises(TypeError):
        Reaction("r", ["A"], [], [1], [], kf=1.0 / unit.sec)


def test_non_species_product_raises():
    """Passing a non-Species object as a product raises TypeError."""
    A = Species("A")
    with pytest.raises(TypeError):
        Reaction("r", [A], ["B"], [1], [1], kf=1.0 / unit.sec)


def test_negative_rate_raises():
    """A negative rate constant is rejected with ValueError."""
    A = Species("A")
    with pytest.raises(ValueError):
        Reaction("r", [A], [], [1], [], kf=-0.1 / unit.sec)


def test_both_rates_zero_warns(caplog):
    """Constructing a reaction with both kf and kr zero emits a warning (but does not raise)."""
    A, B = Species("A"), Species("B")
    with caplog.at_level(logging.WARNING):
        Reaction("r", [A], [B], [1], [1])  # kf=kr=0 by default
    assert any("zero" in rec.message.lower() for rec in caplog.records)


def test_expected_unit_from_order():
    """_expected_unit_from_order returns 1/s, 1/(molar*s), and molar/s for orders 1, 2, and 0."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec)
    assert r._expected_unit_from_order([1]).dimensionality == (1 / unit.sec).dimensionality
    assert r._expected_unit_from_order([1, 1]).dimensionality == (1 / (unit.molar * unit.sec)).dimensionality
    assert r._expected_unit_from_order([]).dimensionality == (unit.molar / unit.sec).dimensionality


def test_validate_rate_scalar_gets_units():
    """_validate_rate assigns the expected unit to a bare scalar and preserves its magnitude."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec)
    expected = r._expected_unit_from_order([1])  # 1/s
    q = r._validate_rate(2.0, expected)
    assert q.dimensionality == expected.dimensionality
    assert q.magnitude == pytest.approx(2.0)


def test_validate_rate_wrong_units_raises():
    """_validate_rate rejects a quantity whose dimensionality doesn't match the expected unit."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec)
    expected = r._expected_unit_from_order([1])  # 1/s
    with pytest.raises(ValueError):
        r._validate_rate(5.0 * unit.molar, expected)  # concentration, not a rate

def test_negative_reverse_rate_raises():
    """A negative reverse rate constant is rejected just like a negative forward one."""
    A, B = Species("A"), Species("B")
    with pytest.raises(ValueError):
        Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec, kr=-0.5 / unit.sec)


def test_species_repr():
    """repr(Species) shows the species ID."""
    assert repr(Species("ATP")) == "Species(ATP)"


def test_display_irreversible_with_stoichiometry():
    """display() writes stoichiometries above 1, joins species with '+', and lists kf."""
    A, B, C, D = Species("A"), Species("B"), Species("C"), Species("D")
    r = Reaction("r", [A, B], [C, D], [1, 2], [3, 1], kf=3.0 / unit.sec)
    assert r.display() == "A + 2 B ---> 3 C + D  // kf = 3.0 / second"


def test_repr_of_reversible_reaction():
    """repr(Reaction) wraps display(); a reversible reaction gets a double arrow and both rates."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [2], kf=1.0 / unit.sec, kr=0.5 / unit.sec)
    assert repr(r) == "Reaction(A <---> 2 B  // kr = 0.5 / second // kf = 1.0 / second)"


def test_display_reverse_only():
    """With kf = 0 only the reverse arrow head and kr are shown."""
    A, B = Species("A"), Species("B")
    s = Reaction("r", [A], [B], [1], [1], kr=0.5 / unit.sec).display()
    assert s.startswith("A <---")
    assert ">" not in s
    assert "kr = 0.5 / second" in s
    assert "kf" not in s


def test_display_zero_order_synthesis():
    """A synthesis reaction with no reactants starts with the arrow."""
    A = Species("A")
    r = Reaction("syn", [], [A], [], [1], kf=2.0 / unit.sec)
    assert r.display() == "---> A  // kf = 2.0 / second"


def test_expected_unit_counts_stoichiometry():
    """Reaction order is the summed stoichiometry: [2] matches [1, 1], and order 3 gives 1/(M^2 s)."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec)
    assert r._expected_unit_from_order([2]) == r._expected_unit_from_order([1, 1])
    third = r._expected_unit_from_order([1, 1, 1])
    assert third.dimensionality == (1 / (unit.molar**2 * unit.sec)).dimensionality


def test_validate_rate_converts_compatible_units():
    """_validate_rate converts a compatible quantity into the expected unit."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec)
    q1 = r._validate_rate(6.0 / unit.minute, r._expected_unit_from_order([1]))
    assert q1.units == unit.second**-1
    assert q1.magnitude == pytest.approx(0.1)
    q2 = r._validate_rate(1e-3 / (unit.millimolar * unit.sec), r._expected_unit_from_order([1, 1]))
    assert q2.units == (1 / (unit.molar * unit.second)).units
    assert q2.magnitude == pytest.approx(1.0)


@pytest.mark.parametrize("rate", [-2.0, -2.0 / unit.sec])
def test_validate_rate_negative_raises(rate):
    """_validate_rate rejects a negative rate, bare or with units."""
    A, B = Species("A"), Species("B")
    r = Reaction("r", [A], [B], [1], [1], kf=1.0 / unit.sec)
    with pytest.raises(ValueError, match="negative"):
        r._validate_rate(rate, r._expected_unit_from_order([1]))
