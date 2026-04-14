"""Tests for ``_snapshot_structure`` — the first-run id snapshot."""

import pytest

import swapasym_cmd as cmd
from _fakes import FakeStructure, residues_from_pairs


def test_snapshot_populates_both_attrs():
    residues = residues_from_pairs([("A", "A"), ("A", "E"), ("B", "B")])
    structure = FakeStructure(residues)

    cmd._snapshot_structure(structure)

    assert residues[0].auth_asym_id == "A"
    assert residues[0].label_asym_id == "A"
    assert residues[1].auth_asym_id == "A"
    assert residues[1].label_asym_id == "E"
    assert residues[2].auth_asym_id == "B"
    assert residues[2].label_asym_id == "B"
    assert getattr(structure, cmd.SNAPSHOT_FLAG) is True


def test_snapshot_is_idempotent_across_chain_id_edits():
    """After a swap, a second snapshot must NOT overwrite the originals."""
    residues = residues_from_pairs([("A", "E"), ("A", "F")])
    structure = FakeStructure(residues)

    cmd._snapshot_structure(structure)
    # Simulate a swap having occurred.
    for r in residues:
        r.chain_id = r.label_asym_id
    # Call snapshot again — it must short-circuit on the flag.
    cmd._snapshot_structure(structure)

    # auth side still the original authors' chain_id, not the swapped one.
    assert [r.auth_asym_id for r in residues] == ["A", "A"]
    assert [r.label_asym_id for r in residues] == ["E", "F"]


def test_snapshot_preserves_empty_label_values():
    """An empty mmcif_chain_id on a subset of residues is still recorded
    verbatim — the caller (``_apply_side``) handles the skip explicitly."""
    residues = residues_from_pairs([("A", "A"), ("A", ""), ("A", "E")])
    structure = FakeStructure(residues)

    cmd._snapshot_structure(structure)

    assert [r.label_asym_id for r in residues] == ["A", "", "E"]


def test_snapshot_raises_not_mmcif_sentinel_when_no_mmcif_data():
    """Structures loaded from .pdb have no mmcif_chain_id on any residue;
    the bundle raises the internal ``_NotMmcif`` sentinel (not a UserError)
    so auto-populate callers can skip silently while explicit callers
    translate it into a user-visible UserError."""
    residues = residues_from_pairs([("A", ""), ("A", ""), ("B", "")])
    structure = FakeStructure(residues)

    with pytest.raises(cmd._NotMmcif):
        cmd._snapshot_structure(structure)


def test_snapshot_accepts_empty_structure():
    """Empty structure shouldn't raise; nothing to snapshot."""
    structure = FakeStructure([])

    cmd._snapshot_structure(structure)

    assert getattr(structure, cmd.SNAPSHOT_FLAG) is True
