"""Test fixtures: stub the chimerax.* imports and load ``src/cmd.py``
under a non-shadowing name (``swapasym_cmd``).

``cmd`` as a bare module name clashes with Python's stdlib ``cmd`` (used
transitively by ``pdb``), so tests import the bundle's command module as
``swapasym_cmd`` via importlib instead of adding ``src/`` to ``sys.path``.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest


_BUNDLE_ROOT = Path(__file__).resolve().parent.parent
_CMD_PATH = _BUNDLE_ROOT.joinpath("src", "cmd.py")


class _StubAtomicStructure:
    """Minimal stand-in for chimerax.atomic.AtomicStructure (isinstance target)."""

    register_attr = MagicMock()


class _StubUserError(Exception):
    """Stand-in for chimerax.core.errors.UserError."""


def _install_stubs() -> None:
    fake_chimerax = MagicMock()
    fake_core = MagicMock()
    fake_atomic = MagicMock()
    fake_commands = MagicMock()
    fake_errors = MagicMock()
    fake_models = MagicMock()
    fake_models.ADD_MODELS = "add models"
    fake_logger = MagicMock()
    fake_logger.html_table_params = "border=1 cellpadding=4 cellspacing=0"

    fake_atomic.Residue = MagicMock()
    fake_atomic.Residue.register_attr = MagicMock()
    fake_atomic.AtomicStructure = _StubAtomicStructure
    fake_atomic.AtomicStructuresArg = MagicMock()

    fake_commands.CmdDesc = MagicMock()
    fake_commands.EnumOf = MagicMock()
    fake_commands.BoolArg = MagicMock()
    fake_commands.run = MagicMock()

    fake_errors.UserError = _StubUserError

    sys.modules["chimerax"] = fake_chimerax
    sys.modules["chimerax.core"] = fake_core
    sys.modules["chimerax.core.commands"] = fake_commands
    sys.modules["chimerax.core.errors"] = fake_errors
    sys.modules["chimerax.core.models"] = fake_models
    sys.modules["chimerax.core.logger"] = fake_logger
    sys.modules["chimerax.atomic"] = fake_atomic


def _load_cmd_module():
    spec = importlib.util.spec_from_file_location("swapasym_cmd", _CMD_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["swapasym_cmd"] = module
    spec.loader.exec_module(module)
    return module


_install_stubs()
_load_cmd_module()


# Make the bundle's test helpers (``_fakes``) importable without installing a package.
_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))


@pytest.fixture(autouse=True)
def reset_module_state():
    import swapasym_cmd as cmd

    cmd._attrs_registered_sessions.clear()
    cmd._add_models_handlers.clear()
    cmd.Residue.register_attr.reset_mock()
    cmd.AtomicStructure.register_attr.reset_mock()
    cmd.run.reset_mock()
    yield
    cmd._attrs_registered_sessions.clear()
    cmd._add_models_handlers.clear()
