"""ChimeraX-SwapAsym bundle entry point."""

from chimerax.core.toolshed import BundleAPI


class _SwapasymAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """Register custom attributes at bundle-init time.

        Attribute registration must happen before any session restore path
        tries to populate residue/structure attributes on this session;
        otherwise session-persistent values are silently dropped.
        """
        from . import cmd

        cmd._register_attrs(session)

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
