"""ChimeraX-SwapAsym bundle entry point."""

from chimerax.core.toolshed import BundleAPI


class _SwapasymAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """Register attributes and subscribe to ADD_MODELS.

        Every structure opened after bundle init gets ``auth_asym_id`` and
        ``label_asym_id`` populated automatically (for mmCIF structures),
        so residue-attribute selectors like ``::label_asym_id="E"`` work
        without having to run ``swapasym`` first.
        """
        from . import cmd

        cmd.install(session)

    @staticmethod
    def finish(session, bundle_info):
        from . import cmd

        cmd.uninstall(session)

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from chimerax.core.commands import register

        from . import cmd

        if command_info.name != "swapasym":
            logger.warning(
                "SwapAsym bundle: unexpected command name "
                f"{command_info.name!r}; not registered. "
                "Check pyproject.toml command declaration."
            )
            return

        register(command_info.name, cmd.swapasym_desc, cmd.swapasym)


bundle_api = _SwapasymAPI()
