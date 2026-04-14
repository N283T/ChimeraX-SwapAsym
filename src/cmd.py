"""Command implementation for ChimeraX-SwapAsym.

Swap between auth_asym_id (PDB chain) and label_asym_id (mmCIF chain) for
chain_id on atomic structures. ``Residue.mmcif_chain_id`` is read-only at
the C binding level, so the bundle snapshots both IDs into custom Residue
attributes (``auth_asym_id`` / ``label_asym_id``) on the first run for each
structure and then reassigns ``Residue.chain_id`` to the requested side.

The custom attributes are registered via ``Residue.register_attr`` so they
become addressable from ChimeraX atom-specs, e.g. ``sel ::label_asym_id="E"``.
"""

import traceback
from weakref import WeakSet

from chimerax.atomic import AtomicStructure, AtomicStructuresArg, Residue
from chimerax.core.commands import CmdDesc, EnumOf
from chimerax.core.errors import UserError

AUTH_ATTR = "auth_asym_id"
LABEL_ATTR = "label_asym_id"
SNAPSHOT_FLAG = "_swapasym_snapshotted"
MODES = ("auto", "label", "auth")

# Sessions on which register_attr has already been invoked. Keyed per-session
# because ChimeraX supports multiple sessions per process; registration is
# also what makes the attributes survive session save/load.
_attrs_registered_sessions: WeakSet = WeakSet()

# Session -> ADD_MODELS handler id, so ``uninstall`` can remove the handler
# when the bundle is unloaded. Regular dict because ChimeraX session objects
# are not weakly referenceable on all platforms; entries are removed in
# ``uninstall``.
_add_models_handlers: dict = {}


def _register_attrs(session):
    """Register custom Residue and Structure attributes once per session.

    Registration is required for session persistence: the custom residue
    attributes and the per-structure snapshot flag survive ``save session``
    / ``open session`` only if they have been registered through ChimeraX's
    attribute registry for the target class.
    """
    if session in _attrs_registered_sessions:
        return
    Residue.register_attr(session, AUTH_ATTR, "SwapAsym", attr_type=str)
    Residue.register_attr(session, LABEL_ATTR, "SwapAsym", attr_type=str)
    AtomicStructure.register_attr(session, SNAPSHOT_FLAG, "SwapAsym", attr_type=bool)
    _attrs_registered_sessions.add(session)


def _try_populate(session, model):
    """Snapshot auth/label ids on a newly-added model, silent on non-mmCIF.

    ``_snapshot_structure`` raises UserError for structures loaded without
    mmcif_chain_id data (plain .pdb, for example). In the auto-populate
    context we want that to be a silent no-op: the bundle should not
    complain about every non-mmCIF file the user opens. When ``swapasym``
    is later invoked explicitly on such a structure, the same UserError
    will surface with its full message.
    """
    if not isinstance(model, AtomicStructure):
        return
    try:
        _snapshot_structure(model)
    except UserError:
        return
    except Exception:
        session.logger.warning(
            f"swapasym: failed to auto-populate {getattr(model, 'name', model)}:\n"
            f"{traceback.format_exc()}"
        )


def _on_add_models(session, trigger_name, models):
    for model in models:
        _try_populate(session, model)


def install(session):
    """Register attributes and subscribe to ADD_MODELS.

    Called from ``BundleAPI.initialize``. After install, every newly
    opened atomic structure gets auth_asym_id / label_asym_id populated
    automatically, so ``sel ::label_asym_id="E"`` and the like work
    without the user having to run ``swapasym`` first.
    """
    _register_attrs(session)
    if session in _add_models_handlers:
        return
    from chimerax.core.models import ADD_MODELS

    handler_id = session.triggers.add_handler(
        ADD_MODELS,
        lambda trigger_name, models: _on_add_models(session, trigger_name, models),
    )
    _add_models_handlers[session] = handler_id

    # Populate any structures that were already open when the bundle
    # initializes (e.g. session restore or manual re-init).
    for model in list(session.models):
        _try_populate(session, model)


def uninstall(session):
    """Remove the ADD_MODELS subscription. Called from BundleAPI.finish."""
    handler_id = _add_models_handlers.pop(session, None)
    if handler_id is None:
        return
    try:
        session.triggers.remove_handler(handler_id)
    except Exception:
        session.logger.info(
            f"swapasym: remove_handler failed during uninstall:\n"
            f"{traceback.format_exc()}"
        )


def _iter_target_structures(session, structures):
    """Return atomic structures to operate on.

    Distinguishes "user spec matched nothing" from "no models open at all"
    so the resulting UserError can point to the actual cause.
    """
    if structures is not None:
        targets = list(structures)
        if not targets:
            raise UserError(
                "swapasym: the given atom spec matched no atomic structures"
            )
        return targets
    atomic = [m for m in session.models if isinstance(m, AtomicStructure)]
    if not atomic:
        raise UserError(
            "swapasym: no atomic structures are open. Open a .cif file first."
        )
    return atomic


def _snapshot_structure(structure):
    """Populate auth/label custom attrs on every residue once per structure.

    Uses a structure-level boolean flag to detect first-run; avoids the
    ambiguity of "unset vs empty string" on the residue attributes themselves.
    Raises UserError if the structure has no mmcif_chain_id info (i.e. was
    not loaded from mmCIF), since without label values ``swapasym mode label``
    would be a no-op.
    """
    if getattr(structure, SNAPSHOT_FLAG, False):
        return

    missing_label = 0
    for residue in structure.residues:
        if not residue.mmcif_chain_id:
            missing_label += 1
    if missing_label == len(structure.residues) and missing_label > 0:
        raise UserError(
            f"swapasym: {structure} has no mmcif_chain_id values "
            "(structure was not loaded from mmCIF). "
            "Reload from a .cif file to use swapasym."
        )

    for residue in structure.residues:
        residue.auth_asym_id = residue.chain_id
        residue.label_asym_id = residue.mmcif_chain_id
    setattr(structure, SNAPSHOT_FLAG, True)


def _current_mode(structure):
    """Detect whether chain_id matches the auth side, label side, or neither.

    Returns one of:
        "empty"      - structure has zero residues
        "identical"  - auth_asym_id == label_asym_id for every residue
                       (swap is a no-op but not an error)
        "auth"       - every residue has chain_id == auth_asym_id
        "label"      - every residue has chain_id == label_asym_id
        "mixed"      - some residues are on each side, or neither
    """
    residues = structure.residues
    if len(residues) == 0:
        return "empty"

    matches_auth = True
    matches_label = True
    all_identical = True
    for residue in residues:
        auth = residue.auth_asym_id
        label = residue.label_asym_id
        cid = residue.chain_id
        if matches_auth and cid != auth:
            matches_auth = False
        if matches_label and cid != label:
            matches_label = False
        if all_identical and auth != label:
            all_identical = False
        if not matches_auth and not matches_label:
            return "mixed"
    if all_identical:
        return "identical"
    if matches_auth:
        return "auth"
    if matches_label:
        return "label"
    return "mixed"


def _apply_side(structure, target_attr):
    """Rewrite chain_id of every residue from the given custom attribute.

    Polymer residues cannot have ``chain_id`` assigned directly at the C
    binding level (ChimeraX raises ``RuntimeError``: "Cannot set polymeric
    chain ID directly from Residue; must use Chain"). For those we batch
    the renames through ``Structure.change_chain_ids``. Non-polymer
    residues (ligands, waters) are written per-residue as before.

    Returns a (changed, skipped) tuple. ``skipped`` covers residues whose
    target attribute is empty (e.g. waters in an mmCIF where the author did
    not assign a label_asym_id); their chain_id is left unchanged and the
    caller is expected to warn the user.
    """
    polymer_targets: dict = {}
    nonpoly_changed = 0
    skipped = 0

    for residue in structure.residues:
        new_cid = getattr(residue, target_attr, None)
        if not new_cid:
            skipped += 1
            continue
        chain = residue.chain
        if chain is None:
            if residue.chain_id != new_cid:
                residue.chain_id = new_cid
                nonpoly_changed += 1
        else:
            if chain.chain_id != new_cid:
                polymer_targets[chain] = new_cid

    polymer_changed = 0
    if polymer_targets:
        chains_list = list(polymer_targets.keys())
        ids_list = [polymer_targets[c] for c in chains_list]
        structure.change_chain_ids(chains_list, ids_list, non_polymeric=False)
        polymer_changed = sum(c.num_residues for c in chains_list)

    return nonpoly_changed + polymer_changed, skipped


def _resolve_target(current, mode):
    """Pick the target side given the detected current state and mode arg.

    ``mode == "auto"`` toggles: auth -> label, label -> auth, identical -> no-op,
    mixed -> UserError (requires explicit mode to avoid silent coercion).
    """
    if mode in ("auth", "label"):
        return mode
    # mode == "auto"
    if current == "auth":
        return "label"
    if current == "label":
        return "auth"
    if current == "identical":
        return "label"
    if current == "empty":
        return "label"
    # current == "mixed"
    raise UserError(
        "swapasym: structure is in a mixed state (chain_id differs from both "
        "auth_asym_id and label_asym_id for some residues). "
        "Re-run with 'mode auth' or 'mode label' to force a side."
    )


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

    for structure in targets:
        _snapshot_structure(structure)
        current = _current_mode(structure)
        num_chains_before = structure.num_chains

        if current == "identical":
            session.logger.warning(
                f"swapasym: {structure} has identical auth_asym_id and "
                "label_asym_id for every residue; swap is a no-op."
            )

        target = _resolve_target(current, mode)
        target_attr = LABEL_ATTR if target == "label" else AUTH_ATTR
        changed, skipped = _apply_side(structure, target_attr)

        if skipped:
            session.logger.warning(
                f"swapasym: {structure} skipped {skipped} residues with "
                f"empty {target_attr} (left on previous side)."
            )

        session.logger.info(
            f"swapasym: {structure} {current} -> {target} "
            f"({changed}/{len(structure.residues)} residues changed, "
            f"{num_chains_before} -> {structure.num_chains} chains)"
        )


swapasym_desc = CmdDesc(
    optional=[("structures", AtomicStructuresArg)],
    keyword=[("mode", EnumOf(MODES))],
    synopsis="Swap chain_id between auth_asym_id and label_asym_id",
)
