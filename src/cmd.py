"""Command implementation for ChimeraX-SwapAsym.

Swap between auth_asym_id (PDB chain) and label_asym_id (mmCIF chain) for
chain_id on atomic structures. Since ``Residue.mmcif_chain_id`` is read-only
at the C binding level, we snapshot both IDs into custom Residue attributes
(``auth_asym_id`` / ``label_asym_id``) on first run, then reassign
``Residue.chain_id`` to the requested side.
"""

from chimerax.atomic import AtomicStructuresArg, Residue
from chimerax.core.commands import CmdDesc, EnumOf
from chimerax.core.errors import UserError

AUTH_ATTR = "auth_asym_id"
LABEL_ATTR = "label_asym_id"
MODES = ("auto", "label", "auth")

_attrs_registered = False


def _register_attrs(session):
    """Register custom Residue attributes once per ChimeraX session."""
    global _attrs_registered
    if _attrs_registered:
        return
    Residue.register_attr(session, AUTH_ATTR, "SwapAsym", attr_type=str)
    Residue.register_attr(session, LABEL_ATTR, "SwapAsym", attr_type=str)
    _attrs_registered = True


def _iter_target_structures(session, structures):
    """Yield atomic structures to operate on."""
    if structures is not None:
        return list(structures)
    from chimerax.atomic import AtomicStructure
    return [m for m in session.models if isinstance(m, AtomicStructure)]


def _snapshot_ids(structure):
    """Save current chain_id as auth and mmcif_chain_id as label.

    Only populates the custom attributes on residues where they are missing,
    so repeated swaps preserve the original values from the first call.
    """
    for residue in structure.residues:
        if getattr(residue, AUTH_ATTR, None) in (None, ""):
            residue.auth_asym_id = residue.chain_id
        if getattr(residue, LABEL_ATTR, None) in (None, ""):
            residue.label_asym_id = residue.mmcif_chain_id


def _current_mode(structure):
    """Detect whether chain_id currently matches the auth or label side.

    Returns ``"auth"`` or ``"label"``. If residues are inconsistent, returns
    ``"mixed"``. When the structure has never been snapshotted, the default
    is treated as ``"auth"`` (ChimeraX's default after mmCIF load).
    """
    residues = structure.residues
    if len(residues) == 0:
        return "auth"
    first = residues[0]
    if getattr(first, AUTH_ATTR, None) in (None, ""):
        return "auth"

    matches_auth = True
    matches_label = True
    for residue in residues:
        if matches_auth and residue.chain_id != residue.auth_asym_id:
            matches_auth = False
        if matches_label and residue.chain_id != residue.label_asym_id:
            matches_label = False
        if not matches_auth and not matches_label:
            return "mixed"
    if matches_label:
        return "label"
    return "auth"


def _apply_side(structure, target_attr):
    """Set chain_id of every residue from the given custom attribute."""
    changed = 0
    for residue in structure.residues:
        new_cid = getattr(residue, target_attr, None)
        if not new_cid:
            continue
        if residue.chain_id != new_cid:
            residue.chain_id = new_cid
            changed += 1
    return changed


def swapasym(session, structures=None, mode="auto"):
    """Swap chain_id between auth_asym_id and label_asym_id.

    Parameters
    ----------
    session : chimerax.core.session.Session
    structures : AtomicStructures or None
        Target structures. Defaults to every open atomic structure.
    mode : {"auto", "label", "auth"}
        ``auto`` toggles between the two sides (default). ``label`` forces
        chain_id to the mmCIF label_asym_id. ``auth`` forces chain_id back
        to the original PDB auth_asym_id.
    """
    _register_attrs(session)

    targets = _iter_target_structures(session, structures)
    if not targets:
        raise UserError("swapasym: no atomic structures to operate on")

    for structure in targets:
        _snapshot_ids(structure)
        current = _current_mode(structure)

        if mode == "auto":
            target = "label" if current != "label" else "auth"
        else:
            target = mode

        target_attr = LABEL_ATTR if target == "label" else AUTH_ATTR
        changed = _apply_side(structure, target_attr)

        session.logger.info(
            f"swapasym: {structure} chain_id -> {target}_asym_id "
            f"({changed} residues changed, {structure.num_chains} chains)"
        )


swapasym_desc = CmdDesc(
    optional=[("structures", AtomicStructuresArg)],
    keyword=[("mode", EnumOf(MODES))],
    synopsis="Swap chain_id between auth_asym_id and label_asym_id",
)
