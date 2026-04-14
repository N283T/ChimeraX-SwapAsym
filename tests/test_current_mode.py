"""Tests for ``_current_mode``."""

import swapasym_cmd as cmd
from _fakes import FakeStructure, residues_from_pairs


def _snapshot(pairs):
    structure = FakeStructure(residues_from_pairs(pairs))
    cmd._snapshot_structure(structure)
    return structure


def test_empty_structure_returns_empty():
    assert cmd._current_mode(FakeStructure([])) == "empty"


def test_auth_state():
    structure = _snapshot([("A", "A"), ("A", "E"), ("B", "B")])
    assert cmd._current_mode(structure) == "auth"


def test_label_state_after_swap():
    structure = _snapshot([("A", "A"), ("A", "E"), ("B", "B")])
    for r in structure.residues:
        r.chain_id = r.label_asym_id
    assert cmd._current_mode(structure) == "label"


def test_identical_ids_detected():
    """If auth == label for every residue, swap is a no-op and we say so."""
    structure = _snapshot([("A", "A"), ("B", "B"), ("C", "C")])
    assert cmd._current_mode(structure) == "identical"


def test_mixed_state_when_some_residues_are_swapped_and_some_not():
    # Use residues where auth != label so "on auth side" vs "on label side"
    # is unambiguous.
    structure = _snapshot([("A", "E"), ("A", "F")])
    # Second residue gets flipped to label, first stays on auth -> mixed.
    structure.residues[1].chain_id = "F"
    assert cmd._current_mode(structure) == "mixed"


def test_mixed_state_when_chain_id_matches_neither():
    structure = _snapshot([("A", "E")])
    structure.residues[0].chain_id = "Z"
    assert cmd._current_mode(structure) == "mixed"
