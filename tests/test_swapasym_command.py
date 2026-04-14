"""Integration-ish tests for the top-level ``swapasym`` function.

These exercise the full flow (register -> iterate targets -> snapshot ->
mode detect -> apply -> log) using fake Session/Structure objects.
"""

import pytest

import swapasym_cmd as cmd
from _fakes import FakeSession, FakeStructure, residues_from_pairs


def _session_with(pairs, name="4hhb"):
    structure = FakeStructure(residues_from_pairs(pairs), name=name)
    return FakeSession(models=[structure]), structure


def test_registers_attrs_once_per_session():
    session, _ = _session_with([("A", "E")])

    cmd.swapasym(session)
    cmd.swapasym(session)
    cmd.swapasym(session)

    # Two Residue attrs registered exactly once (idempotent per session).
    assert cmd.Residue.register_attr.call_count == 2
    # One AtomicStructure attr (the snapshot flag).
    assert cmd.AtomicStructure.register_attr.call_count == 1


def test_registers_attrs_again_on_new_session():
    """Custom-attr registration is session-scoped, not process-scoped."""
    session1, _ = _session_with([("A", "E")], name="first")
    session2, _ = _session_with([("A", "E")], name="second")

    cmd.swapasym(session1)
    cmd.swapasym(session2)

    # 2 Residue attrs * 2 sessions.
    assert cmd.Residue.register_attr.call_count == 4
    # 1 AtomicStructure attr * 2 sessions.
    assert cmd.AtomicStructure.register_attr.call_count == 2


def test_toggle_round_trip():
    """swapasym auto twice returns chain_id to the starting state."""
    session, structure = _session_with([("A", "A"), ("A", "E"), ("B", "B")])

    cmd.swapasym(session)
    after_first = [r.chain_id for r in structure.residues]
    cmd.swapasym(session)
    after_second = [r.chain_id for r in structure.residues]

    assert after_first == ["A", "E", "B"]
    assert after_second == ["A", "A", "B"]


def test_explicit_mode_label():
    session, structure = _session_with([("A", "A"), ("A", "E")])

    cmd.swapasym(session, mode="label")

    assert [r.chain_id for r in structure.residues] == ["A", "E"]


def test_explicit_mode_auth_restores_even_from_mixed_state():
    """Mixed state + explicit mode auth -> always restore cleanly."""
    session, structure = _session_with([("A", "A"), ("A", "E")])
    cmd.swapasym(session, mode="label")
    structure.residues[0].chain_id = "Z"  # inject mixed state

    cmd.swapasym(session, mode="auth")

    assert [r.chain_id for r in structure.residues] == ["A", "A"]


def test_mixed_state_with_auto_raises_usererror():
    session, structure = _session_with([("A", "A"), ("A", "E")])
    cmd.swapasym(session, mode="label")
    structure.residues[0].chain_id = "Z"

    with pytest.raises(cmd.UserError):
        cmd.swapasym(session)  # mode=auto default


def test_raises_when_no_atomic_structures_open():
    session = FakeSession(models=[])

    with pytest.raises(cmd.UserError) as exc_info:
        cmd.swapasym(session)
    assert "no atomic structures" in str(exc_info.value).lower()


def test_raises_when_user_spec_matched_nothing():
    session, _ = _session_with([("A", "A")])

    with pytest.raises(cmd.UserError) as exc_info:
        cmd.swapasym(session, structures=[])
    assert "matched no atomic structures" in str(exc_info.value).lower()


def test_structure_loaded_from_pdb_raises_on_first_use():
    """No mmcif_chain_id on any residue -> UserError (can't swap to label)."""
    session, _ = _session_with([("A", ""), ("B", "")], name="pdb_only")

    with pytest.raises(cmd.UserError) as exc_info:
        cmd.swapasym(session)
    assert "mmcif" in str(exc_info.value).lower()


def test_warns_when_identical_ids():
    """Structure whose auth == label for every residue: log a warning so
    the user understands why swap appears to do nothing."""
    session, _ = _session_with([("A", "A"), ("B", "B")])

    cmd.swapasym(session)

    warnings = " ".join(session.logger.warning_msgs).lower()
    assert "identical" in warnings


def test_warns_when_residues_skipped_due_to_empty_label():
    session, structure = _session_with([("A", "A"), ("A", ""), ("A", "E")])

    cmd.swapasym(session, mode="label")

    warnings = " ".join(session.logger.warning_msgs).lower()
    assert "skipped" in warnings
    # One residue has empty label, so one is skipped.
    assert "1 residues" in warnings or "1 residue" in warnings


def test_info_log_is_html_table_with_summary_rows():
    """Success log must be emitted as HTML (is_html=True) and include
    rows for residues / polymer chains / unique chain_ids."""
    session, _ = _session_with([("A", "A"), ("A", "E"), ("A", "F"), ("B", "B")])

    cmd.swapasym(session, mode="label")

    assert len(session.logger.html_info_msgs) == 1
    html = session.logger.html_info_msgs[0]
    assert "<table" in html
    assert "auth &rarr; label" in html
    assert "residues" in html
    assert "polymer chains" in html
    assert "unique chain_ids" in html


def test_info_log_added_table_has_type_column_and_cxcmd_link():
    """When swapping to label, newly-appearing chain_ids should appear in
    a second HTML table row tagged 'added', with a residue-type cell and
    a clickable ``cxcmd:select`` link for the chain id."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True, name="VAL"),
        FakeResidue("A", "A", polymer=True, name="ALA"),
        FakeResidue("A", "E", polymer=False, name="HEM"),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")

    html = session.logger.html_info_msgs[0]
    assert "added" in html
    assert ">HEM<" in html  # residue-type cell
    assert 'href="cxcmd:select #1/E"' in html


def test_info_log_removed_row_when_swapping_back_to_auth():
    """Reverse swap: the label-only ids should appear in a 'removed' row
    (without a cxcmd link, since the id no longer exists post-swap)."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True, name="VAL"),
        FakeResidue("A", "E", polymer=False, name="HEM"),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")  # auth -> label
    session.logger.html_info_msgs.clear()
    session.logger.info_msgs.clear()
    cmd.swapasym(session, mode="auth")  # label -> auth

    html = session.logger.html_info_msgs[0]
    assert "removed" in html
    assert ">HEM<" in html
    # Removed chain_ids are no longer present on the current side, so no link.
    assert 'href="cxcmd:select #1/E"' not in html


def test_multi_structure_all_processed():
    s1 = FakeStructure(residues_from_pairs([("A", "A"), ("A", "E")]), name="s1")
    s2 = FakeStructure(residues_from_pairs([("B", "B"), ("B", "F")]), name="s2")
    session = FakeSession(models=[s1, s2])

    cmd.swapasym(session, mode="label")

    assert [r.chain_id for r in s1.residues] == ["A", "E"]
    assert [r.chain_id for r in s2.residues] == ["B", "F"]


def test_color_option_runs_color_bychain_on_each_target():
    """With color=True, swapasym runs `color <spec> bychain` covering every
    affected structure after the swap."""
    s1 = FakeStructure(residues_from_pairs([("A", "E")]), name="s1", atomspec="#1")
    s2 = FakeStructure(residues_from_pairs([("B", "F")]), name="s2", atomspec="#2")
    session = FakeSession(models=[s1, s2])

    cmd.swapasym(session, mode="label", color=True)

    assert cmd.run.call_count == 1
    args, kwargs = cmd.run.call_args
    # positional: (session, command_text); kwarg: log=False
    assert args[0] is session
    assert args[1] == "color #1 #2 bychain"
    assert kwargs.get("log") is False


def test_color_option_defaults_to_off():
    """Existing invocations without ``color`` must NOT re-color the scene."""
    session, _ = _session_with([("A", "E")])

    cmd.swapasym(session, mode="label")

    cmd.run.assert_not_called()


def test_color_option_skipped_when_no_targets_match():
    """No ``color`` command is issued if the swap itself raised (no targets)."""
    session = FakeSession(models=[])
    with pytest.raises(cmd.UserError):
        cmd.swapasym(session, color=True)
    cmd.run.assert_not_called()
