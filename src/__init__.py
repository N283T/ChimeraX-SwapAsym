"""ChimeraX-SwapAsym - ChimeraX SwapAsym bundle"""

from chimerax.core.toolshed import BundleAPI


class _SwapasymAPI(BundleAPI):
    """Bundle API for ChimeraX-SwapAsym."""

    api_version = 1

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        """Register the 'swapasym' command."""
        from . import cmd

        if command_info.name == "swapasym":
            from chimerax.core.commands import register

            register(
                command_info.name,
                cmd.swapasym_desc,
                cmd.swapasym,
            )


bundle_api = _SwapasymAPI()
