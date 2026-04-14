"""Lightweight fake Residue / Structure / Session used across tests."""

from __future__ import annotations

from typing import Optional

# ``chimerax.atomic.AtomicStructure`` is already stubbed by conftest by the
# time this module is first imported. Using it as a base class means the
# ``isinstance(m, AtomicStructure)`` filter in ``_iter_target_structures``
# accepts instances of ``FakeStructure``.
from chimerax.atomic import AtomicStructure as _AtomicStructure


class FakeResidue:
    """Stand-in for chimerax.atomic.Residue.

    Models auth side (``chain_id``), label side (``mmcif_chain_id``),
    and the bundle-managed custom attributes (``auth_asym_id``,
    ``label_asym_id``) the same way real ChimeraX objects do.

    The ``chain`` back-reference is set by ``FakeStructure`` for polymer
    residues (mirrors ``chimerax.atomic.Residue.chain``); non-polymer
    residues leave it ``None``.
    """

    def __init__(
        self,
        chain_id: str,
        mmcif_chain_id: str = "",
        *,
        polymer: bool = False,
        name: str = "UNK",
    ):
        self.chain_id = chain_id
        self.mmcif_chain_id = mmcif_chain_id
        self.name = name
        self._polymer = polymer
        self.chain: Optional["FakeChain"] = None

    def __repr__(self) -> str:  # pragma: no cover - debug aid only
        return (
            f"FakeResidue(chain_id={self.chain_id!r}, "
            f"mmcif_chain_id={self.mmcif_chain_id!r}, "
            f"polymer={self._polymer!r})"
        )


class FakeChain:
    """Stand-in for chimerax.atomic.Chain.

    Tracks its member polymer residues so ``Structure.change_chain_ids``
    can push a new ``chain_id`` down to each of them in a single shot.
    """

    def __init__(self, chain_id: str, residues: list[FakeResidue]):
        self.chain_id = chain_id
        self._residues = list(residues)
        for residue in self._residues:
            residue.chain = self

    @property
    def num_residues(self) -> int:
        return len(self._residues)

    def _rename(self, new_cid: str) -> None:
        self.chain_id = new_cid
        for residue in self._residues:
            # Bypass the polymer restriction in the fake by assigning the
            # attribute directly — the real C layer updates these too.
            residue.chain_id = new_cid


class FakeStructure(_AtomicStructure):
    """Stand-in for chimerax.atomic.AtomicStructure.

    Inherits from the stubbed ``AtomicStructure`` so ``isinstance`` checks in
    the bundle code succeed. ``num_chains`` is recomputed from the unique
    residue chain_ids so the value reflects the current state.
    """

    def __init__(
        self,
        residues: list[FakeResidue],
        name: str = "fake",
        atomspec: str = "#1",
    ):
        self.residues = list(residues)
        self.name = name
        self.atomspec = atomspec
        # Build one FakeChain per distinct chain_id for all polymer residues.
        polymer_groups: dict[str, list[FakeResidue]] = {}
        for residue in residues:
            if residue._polymer:
                polymer_groups.setdefault(residue.chain_id, []).append(residue)
        self._chains = [
            FakeChain(cid, members) for cid, members in polymer_groups.items()
        ]

    @property
    def num_chains(self) -> int:
        return len(self._chains)

    @property
    def chains(self) -> list[FakeChain]:
        return list(self._chains)

    def change_chain_ids(self, chains, chain_ids, non_polymeric: bool = True):
        for chain, new_cid in zip(chains, chain_ids):
            chain._rename(new_cid)

    def __repr__(self) -> str:  # pragma: no cover - debug aid only
        return f"FakeStructure({self.name!r}, {len(self.residues)} residues)"


class FakeLogger:
    def __init__(self):
        self.info_msgs: list[str] = []
        self.warning_msgs: list[str] = []

    def info(self, msg: str) -> None:
        self.info_msgs.append(msg)

    def warning(self, msg: str) -> None:
        self.warning_msgs.append(msg)


class FakeTriggers:
    """Minimal stand-in for ``session.triggers`` used by install/uninstall.

    Deliberately does NOT catch handler exceptions in ``fire`` (unlike real
    ChimeraX's ``Triggers.activate_trigger``): propagating exceptions in
    tests surfaces real bugs. Set ``remove_raises`` to make the next
    ``remove_handler`` call raise the given exception — used to exercise
    the uninstall failure path.
    """

    def __init__(self):
        self._handlers: dict[int, tuple] = {}
        self._next_id = 1
        self.remove_raises: Exception | None = None

    def add_handler(self, trigger_name: str, callback):
        hid = self._next_id
        self._next_id += 1
        self._handlers[hid] = (trigger_name, callback)
        return hid

    def remove_handler(self, handler_id) -> None:
        if self.remove_raises is not None:
            raise self.remove_raises
        self._handlers.pop(handler_id, None)

    def fire(self, trigger_name: str, payload):
        for tname, callback in list(self._handlers.values()):
            if tname == trigger_name:
                callback(tname, payload)


class FakeSession:
    def __init__(self, models: Optional[list[FakeStructure]] = None):
        self.models = list(models or [])
        self.logger = FakeLogger()
        self.triggers = FakeTriggers()


def residues_from_pairs(
    pairs: list[tuple[str, str]], *, polymer: bool = False
) -> list[FakeResidue]:
    """Shorthand: build a residue list from (chain_id, mmcif_chain_id) pairs.

    When ``polymer=True``, every residue is marked as polymer so the bundle
    code routes ``chain_id`` writes through ``change_chain_ids`` for them.
    """
    return [FakeResidue(cid, lbl, polymer=polymer) for cid, lbl in pairs]
