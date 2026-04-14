"""Tests for ``_resolve_target`` — the toggle decision logic."""

import pytest

import swapasym_cmd as cmd


def test_auto_auth_to_label():
    assert cmd._resolve_target("auth", "auto") == "label"


def test_auto_label_to_auth():
    assert cmd._resolve_target("label", "auto") == "auth"


def test_auto_identical_defaults_to_label():
    # identical is a no-op either way; picking one deterministically.
    assert cmd._resolve_target("identical", "auto") == "label"


def test_auto_empty_defaults_to_label():
    assert cmd._resolve_target("empty", "auto") == "label"


def test_auto_mixed_raises():
    """Mixed state must not silently coerce to either side."""
    with pytest.raises(cmd.UserError) as exc_info:
        cmd._resolve_target("mixed", "auto")
    assert "mixed state" in str(exc_info.value).lower()


@pytest.mark.parametrize("current", ["auth", "label", "mixed", "identical", "empty"])
def test_explicit_label_always_returns_label(current):
    assert cmd._resolve_target(current, "label") == "label"


@pytest.mark.parametrize("current", ["auth", "label", "mixed", "identical", "empty"])
def test_explicit_auth_always_returns_auth(current):
    assert cmd._resolve_target(current, "auth") == "auth"
