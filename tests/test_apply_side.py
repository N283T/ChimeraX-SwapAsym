"""Tests for ``_apply_side`` — the per-residue chain_id rewrite."""

import swapasym_cmd as cmd
from _fakes import FakeStructure, residues_from_pairs


def _snapshot(pairs):
    structure = FakeStructure(residues_from_pairs(pairs))
    cmd._snapshot_structure(structure)
    return structure


def test_apply_label_rewrites_chain_id():
    structure = _snapshot([("A", "A"), ("A", "E"), ("B", "B")])

    changed, skipped = cmd._apply_side(structure, cmd.LABEL_ATTR)

    assert changed == 1  # only residue[1] actually changes (A -> E)
    assert skipped == 0
    assert [r.chain_id for r in structure.residues] == ["A", "E", "B"]


def test_apply_auth_restores_original():
    structure = _snapshot([("A", "A"), ("A", "E"), ("B", "B")])
    # First swap to label.
    cmd._apply_side(structure, cmd.LABEL_ATTR)
    # Then back.
    changed, skipped = cmd._apply_side(structure, cmd.AUTH_ATTR)

    assert skipped == 0
    assert [r.chain_id for r in structure.residues] == ["A", "A", "B"]


def test_apply_side_reports_skipped_when_label_missing():
    """A residue with empty label_asym_id must be counted as skipped,
    not silently left on auth while pretending it succeeded."""
    structure = _snapshot([("A", "A"), ("A", ""), ("A", "E")])

    changed, skipped = cmd._apply_side(structure, cmd.LABEL_ATTR)

    assert skipped == 1
    assert changed == 1  # only the one with a real label
    # Residue with empty label retains its original chain_id.
    assert structure.residues[1].chain_id == "A"


def test_apply_side_is_stable_under_repeat():
    """Re-running the same side is a no-op."""
    structure = _snapshot([("A", "E")])
    cmd._apply_side(structure, cmd.LABEL_ATTR)

    changed, skipped = cmd._apply_side(structure, cmd.LABEL_ATTR)

    assert changed == 0
    assert skipped == 0


def test_apply_side_uses_change_chain_ids_for_polymer_residues():
    """Polymer residues cannot have chain_id assigned directly; the bundle
    must route those writes through Structure.change_chain_ids."""
    from _fakes import FakeStructure, residues_from_pairs

    residues = residues_from_pairs([("A", "X"), ("A", "X")], polymer=True)
    structure = FakeStructure(residues)
    cmd._snapshot_structure(structure)

    changed, skipped = cmd._apply_side(structure, cmd.LABEL_ATTR)

    assert skipped == 0
    assert changed == 2
    # After change_chain_ids: every polymer residue carries the new chain id.
    assert [r.chain_id for r in residues] == ["X", "X"]
    # The Chain object's own chain_id reflects the rename as well.
    assert structure.chains[0].chain_id == "X"


def test_apply_side_mixed_polymer_and_nonpolymer():
    """One structure with both polymer and non-polymer residues should
    handle each correctly (change_chain_ids for polymer, direct for others)."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "X", polymer=True),
        FakeResidue("A", "X", polymer=True),
        FakeResidue("A", "E", polymer=False),  # ligand
        FakeResidue("A", "K", polymer=False),  # water
    ]
    structure = FakeStructure(residues)
    cmd._snapshot_structure(structure)

    changed, _ = cmd._apply_side(structure, cmd.LABEL_ATTR)

    assert [r.chain_id for r in residues] == ["X", "X", "E", "K"]
    assert changed == 4
