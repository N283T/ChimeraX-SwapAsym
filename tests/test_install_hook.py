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


# --- Regression tests: failure modes surfaced by the PR review ---


def test_uninstall_remove_handler_failure_keeps_entry_and_warns(monkeypatch):
    """If ``remove_handler`` raises, the session must remain in the handler
    map so that a subsequent ``install`` does not register a duplicate;
    the failure is surfaced at warning level."""
    session = FakeSession()
    cmd.install(session)
    session.triggers.remove_raises = RuntimeError("stale handler")

    cmd.uninstall(session)  # must not raise

    # Entry retained so re-install does not double-subscribe.
    assert session in cmd._add_models_handlers
    assert any("remove_handler failed" in m for m in session.logger.warning_msgs), (
        session.logger.warning_msgs
    )
    assert session.logger.info_msgs == []


def test_reinstall_after_failed_uninstall_does_not_duplicate_handler(monkeypatch):
    """Keep the dict entry on failure so ``install`` is a no-op on the
    already-registered handler instead of adding a second one."""
    session = FakeSession()
    cmd.install(session)
    session.triggers.remove_raises = RuntimeError("stale handler")
    cmd.uninstall(session)

    # Clear the remove_raises so a later (successful) uninstall can clean up.
    session.triggers.remove_raises = None
    cmd.install(session)

    # Only one handler should be live.
    assert len(session.triggers._handlers) == 1


def test_on_add_models_continues_past_failing_model(monkeypatch):
    """A single model that blows up in ``_snapshot_structure`` must not
    abort the rest of the ADD_MODELS batch."""
    session = FakeSession()
    cmd.install(session)
    good = _structure([("A", "E")], name="good")
    bad = _structure([("B", "F")], name="bad")

    real_snapshot = cmd._snapshot_structure

    def _selective_snapshot(structure):
        if structure.name == "bad":
            raise RuntimeError("boom")
        real_snapshot(structure)

    monkeypatch.setattr(cmd, "_snapshot_structure", _selective_snapshot)

    session.triggers.fire("add models", [bad, good])  # must not raise

    # Good model still got populated after bad model failed.
    assert getattr(good, cmd.SNAPSHOT_FLAG, False) is True
    # Warning emitted naming the bad model.
    assert any("bad" in m for m in session.logger.warning_msgs), (
        session.logger.warning_msgs
    )


def test_add_models_mixed_mmcif_and_non_mmcif_populates_only_mmcif(monkeypatch):
    """A batch containing a PDB-only structure and an mmCIF structure must
    populate the mmCIF one; the PDB-only one is silently skipped."""
    session = FakeSession()
    cmd.install(session)
    mmcif = _structure([("A", "E"), ("B", "F")], name="mmcif")
    pdb_only = _structure([("A", ""), ("B", "")], name="pdb")

    session.triggers.fire("add models", [pdb_only, mmcif])

    assert getattr(mmcif, cmd.SNAPSHOT_FLAG, False) is True
    assert getattr(pdb_only, cmd.SNAPSHOT_FLAG, False) is False
    assert session.logger.warning_msgs == []


def test_try_populate_unexpected_usererror_is_reported(monkeypatch):
    """UserErrors that are NOT ``_NotMmcif`` should no longer be swallowed
    silently: they surface as a warning so auto-populate bugs do not hide.
    """
    session = FakeSession()
    cmd.install(session)
    structure = _structure([("A", "E")], name="weird")

    def _raise_unrelated(_):
        raise cmd.UserError("unrelated failure")

    monkeypatch.setattr(cmd, "_snapshot_structure", _raise_unrelated)

    session.triggers.fire("add models", [structure])

    assert any(
        "weird" in m and "auto-populate" in m for m in session.logger.warning_msgs
    ), session.logger.warning_msgs


def test_install_uninstall_are_per_session():
    """Each session keeps its own handler; uninstalling one leaves the
    other's subscription intact."""
    s1 = FakeSession()
    s2 = FakeSession()
    cmd.install(s1)
    cmd.install(s2)

    assert s1 in cmd._add_models_handlers
    assert s2 in cmd._add_models_handlers
    # Each session has its own triggers object; both get exactly one handler.
    assert s1.triggers is not s2.triggers
    assert len(s1.triggers._handlers) == 1
    assert len(s2.triggers._handlers) == 1

    s1_only = _structure([("A", "E")], name="s1_only")
    s1.triggers.fire("add models", [s1_only])
    assert getattr(s1_only, cmd.SNAPSHOT_FLAG, False) is True

    cmd.uninstall(s1)
    assert s1 not in cmd._add_models_handlers
    assert s2 in cmd._add_models_handlers
    assert len(s2.triggers._handlers) == 1


def test_install_uninstall_install_cycle_leaves_single_handler():
    """Plugin reload (install → uninstall → install) must end up with
    exactly one active handler."""
    session = FakeSession()
    cmd.install(session)
    cmd.uninstall(session)
    cmd.install(session)

    assert len(session.triggers._handlers) == 1
    assert session in cmd._add_models_handlers


def test_install_prepopulate_survives_failing_model(monkeypatch):
    """A failing pre-existing model at install time must not abort
    pre-population of the rest."""
    good = _structure([("A", "E")], name="good")
    bad = _structure([("B", "F")], name="bad")
    session = FakeSession(models=[bad, good])

    real_snapshot = cmd._snapshot_structure

    def _selective_snapshot(structure):
        if structure.name == "bad":
            raise RuntimeError("init boom")
        real_snapshot(structure)

    monkeypatch.setattr(cmd, "_snapshot_structure", _selective_snapshot)

    cmd.install(session)  # must not raise

    assert getattr(good, cmd.SNAPSHOT_FLAG, False) is True
    assert any("bad" in m for m in session.logger.warning_msgs)


def test_non_atomic_structure_model_in_payload_is_ignored():
    """ADD_MODELS fires for volumes, surfaces, etc. Non-AtomicStructure
    payload items are silently ignored with no log noise."""

    class _FakeVolume:
        name = "map.mrc"

    session = FakeSession()
    cmd.install(session)

    session.triggers.fire("add models", [_FakeVolume()])  # must not raise
    assert session.logger.warning_msgs == []
