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


def _tfoot_region(html: str) -> str:
    """Extract the ``<tfoot>…</tfoot>`` substring so tests can scope
    assertions to the footer and not accidentally match cells from
    the mapping body."""
    start = html.index("<tfoot>")
    end = html.index("</tfoot>") + len("</tfoot>")
    return html[start:end]


def test_info_log_is_single_html_table_with_tfoot_unique_counts():
    """Success log emits exactly one HTML table; the swap header is in
    ``<thead>``, the groupby mapping in ``<tbody>``, and the per-side
    unique-chain_id counts in ``<tfoot>``."""
    session, _ = _session_with([("A", "A"), ("A", "E"), ("A", "F"), ("B", "B")])

    cmd.swapasym(session, mode="label")

    assert len(session.logger.html_info_msgs) == 1
    html = session.logger.html_info_msgs[0]
    assert html.count("<table") == 1
    # Table skeleton all present; summary header spans all three columns.
    assert "<thead>" in html
    assert "<tbody>" in html
    assert "<tfoot>" in html
    assert 'colspan="3"' in html
    # Header carries the direction text and column labels.
    assert "auth &rarr; label" in html
    assert "auth_asym_id" in html and "label_asym_id" in html
    # Thead column row has exactly three <th> cells (blank middle column),
    # plus the summary <th colspan="3"> on the preceding row.
    thead_region = html[html.index("<thead>") : html.index("</thead>")]
    assert thead_region.count("<th>") == 3  # three unadorned <th> cells
    assert thead_region.count('<th colspan="3">') == 1
    # Tfoot asserts — scope to the footer region so we do not false-match
    # bold cells inside the mapping body.
    tfoot = _tfoot_region(html)
    assert "<i>unique</i>" in tfoot
    # auth_unique (2) appears before label_unique (4) in the footer row.
    assert tfoot.index("<b>2</b>") < tfoot.index("<b>4</b>")
    # Noise from earlier iterations is gone.
    assert "residues changed" not in html
    assert "polymer chains" not in html


def test_info_log_mapping_columns_are_fixed_auth_label():
    """Column order is always ``auth_asym_id | label_asym_id`` regardless
    of swap direction; swap direction is conveyed by the summary header."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True),
        FakeResidue("A", "E", polymer=False),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")
    forward_html = session.logger.html_info_msgs[0]

    session.logger.html_info_msgs.clear()
    session.logger.info_msgs.clear()
    cmd.swapasym(session, mode="auth")
    reverse_html = session.logger.html_info_msgs[0]

    for html in (forward_html, reverse_html):
        assert html.index("auth_asym_id") < html.index("label_asym_id")


def test_info_log_mapping_uses_rowspan_groupby():
    """Both the auth cell and the arrow cell are rowspan-merged to cover
    every label row — so a regression dropping the arrow's rowspan would
    still produce one rowspan=3, but break the groupby visually."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True),
        FakeResidue("A", "E", polymer=False),
        FakeResidue("A", "F", polymer=False),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")

    html = session.logger.html_info_msgs[0]
    # Two rowspan="3" cells (auth + arrow) covering labels A, E, F.
    assert html.count('rowspan="3"') == 2


def test_info_log_clickable_cells_follow_current_chain_id_side_forward():
    """After auth→label, label cells link via ``#1/X`` (current) and auth
    cells link via ``::auth_asym_id='X'`` (other side)."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True),
        FakeResidue("A", "E", polymer=False),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")

    html = session.logger.html_info_msgs[0]
    assert "cxcmd:select ::auth_asym_id=&#39;A&#39;" in html
    assert "cxcmd:select #1/E" in html
    # Auth cell is wrapped in <b>; it must use the custom-attr selector
    # rather than the standard spec after auth→label.
    assert '<b><a href="cxcmd:select #1/' not in html


def test_info_log_clickable_cells_follow_current_chain_id_side_reverse():
    """After label→auth, auth cells link via ``#1/X`` (current) and label
    cells link via ``::label_asym_id='X'`` (other side). Mirror of the
    forward case."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True),
        FakeResidue("A", "E", polymer=False),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")
    session.logger.html_info_msgs.clear()
    session.logger.info_msgs.clear()
    cmd.swapasym(session, mode="auth")

    html = session.logger.html_info_msgs[0]
    assert "cxcmd:select ::label_asym_id=&#39;E&#39;" in html
    # The auth cell (bolded) now uses the standard spec.
    assert '<b><a href="cxcmd:select #1/A">A</a></b>' in html
    # Label side must NOT use the standard spec after label→auth.
    assert "cxcmd:select #1/E" not in html


def test_info_log_clickable_cells_fall_back_to_attr_when_no_atomspec():
    """Structures with empty ``atomspec`` have no valid standard selector;
    both sides must still be clickable via the custom-attr selector."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True),
        FakeResidue("A", "E", polymer=False),
    ]
    structure = FakeStructure(residues, atomspec="")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")

    html = session.logger.html_info_msgs[0]
    # Both sides use custom-attr selectors when atomspec is absent.
    assert "cxcmd:select ::auth_asym_id=&#39;A&#39;" in html
    assert "cxcmd:select ::label_asym_id=&#39;E&#39;" in html
    # The broken ``#/X`` spec must not appear.
    assert "cxcmd:select /A" not in html
    assert "cxcmd:select /E" not in html


def test_info_log_mapping_has_direction_arrow():
    """The mapping table's middle column carries ``&rarr;`` after an
    auth→label swap and ``&larr;`` after a label→auth swap. The centered
    ``text-align:center`` cell style pins the match to the mapping table
    rather than the summary header (which also uses ``&rarr;``)."""
    from _fakes import FakeResidue, FakeStructure

    residues = [
        FakeResidue("A", "A", polymer=True),
        FakeResidue("A", "E", polymer=False),
    ]
    structure = FakeStructure(residues, atomspec="#1")
    session = FakeSession(models=[structure])

    cmd.swapasym(session, mode="label")
    forward = session.logger.html_info_msgs[0]
    assert "text-align:center" in forward
    assert 'text-align:center">&rarr;' in forward
    assert 'text-align:center">&larr;' not in forward

    session.logger.html_info_msgs.clear()
    session.logger.info_msgs.clear()
    cmd.swapasym(session, mode="auth")
    reverse = session.logger.html_info_msgs[0]
    assert 'text-align:center">&larr;' in reverse
    assert 'text-align:center">&rarr;' not in reverse


def test_info_log_has_no_added_removed_tables():
    """After simplification, only the summary + mapping tables remain."""
    session, _ = _session_with([("A", "A"), ("A", "E")])

    cmd.swapasym(session, mode="label")

    html = session.logger.html_info_msgs[0]
    assert "added" not in html
    assert "removed" not in html
    assert "details" not in html


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
