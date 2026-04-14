"""Tests for ``install`` / ``uninstall`` and the ADD_MODELS auto-populate hook."""

import swapasym_cmd as cmd
from _fakes import FakeSession, FakeStructure, residues_from_pairs


def _structure(pairs, name="fake"):
    return FakeStructure(residues_from_pairs(pairs), name=name)


def test_install_subscribes_once_and_registers_attrs():
    session = FakeSession()

    cmd.install(session)
    cmd.install(session)
    cmd.install(session)

    # Handler table gets exactly one entry, registered against ADD_MODELS.
    assert session in cmd._add_models_handlers
    assert len(session.triggers._handlers) == 1
    trigger_name, _ = next(iter(session.triggers._handlers.values()))
    assert trigger_name == "add models"


def test_install_populates_already_open_structures():
    """Structures open before install are snapshotted at install time."""
    structure = _structure([("A", "A"), ("A", "E"), ("B", "B")])
    session = FakeSession(models=[structure])

    cmd.install(session)

    assert getattr(structure, cmd.SNAPSHOT_FLAG, False) is True
    assert structure.residues[0].auth_asym_id == "A"
    assert structure.residues[1].label_asym_id == "E"


def test_newly_added_model_triggers_populate():
    session = FakeSession()
    cmd.install(session)

    structure = _structure([("A", "E"), ("B", "F")])
    session.models.append(structure)
    session.triggers.fire("add models", [structure])

    assert getattr(structure, cmd.SNAPSHOT_FLAG, False) is True
    assert structure.residues[0].label_asym_id == "E"
    assert structure.residues[1].auth_asym_id == "B"


def test_non_mmcif_structure_is_silently_skipped():
    """PDB-only structures have empty mmcif_chain_id on every residue and
    would raise in ``_snapshot_structure``; the hook must swallow it."""
    session = FakeSession()
    cmd.install(session)
    pdb_only = _structure([("A", ""), ("A", ""), ("B", "")], name="pdb_only")

    session.triggers.fire("add models", [pdb_only])  # must not raise

    assert getattr(pdb_only, cmd.SNAPSHOT_FLAG, False) is False
    assert session.logger.warning_msgs == []


def test_uninstall_removes_handler():
    session = FakeSession()
    cmd.install(session)

    cmd.uninstall(session)

    assert session not in cmd._add_models_handlers
    assert session.triggers._handlers == {}


def test_uninstall_is_idempotent():
    session = FakeSession()
    cmd.install(session)
    cmd.uninstall(session)
    cmd.uninstall(session)  # must not raise


def test_uninstall_without_install_is_noop():
    session = FakeSession()
    cmd.uninstall(session)  # must not raise
    assert session.triggers._handlers == {}


def test_snapshot_flag_short_circuits_swapasym_snapshot():
    """After the hook pre-snapshots a structure, ``swapasym`` picks up the
    existing snapshot and does not overwrite auth values."""
    structure = _structure([("A", "E"), ("B", "F")])
    session = FakeSession(models=[structure])
    cmd.install(session)

    # Simulate a later chain_id edit (e.g. the user renamed a chain) and
    # make sure a subsequent snapshot_structure call does not re-run.
    structure.residues[0].chain_id = "Z"
    cmd._snapshot_structure(structure)

    # auth_asym_id was captured at install time, not after the edit.
    assert structure.residues[0].auth_asym_id == "A"
