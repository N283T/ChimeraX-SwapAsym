# ChimeraX-SwapAsym

Swap `chain_id` between `auth_asym_id` (PDB) and `label_asym_id` (mmCIF) on
loaded atomic structures.

## Motivation

mmCIF files carry two chain identifiers:

- `auth_asym_id` — the author / PDB chain ID (what ChimeraX exposes as
  `Residue.chain_id` and uses everywhere in the UI).
- `label_asym_id` — the mmCIF identifier, where each polymer / branched
  entity / ligand / waters group typically gets its own ID (read-only
  `Residue.mmcif_chain_id`).

This bundle copies both into custom Residue attributes on first use so you
can switch the primary `chain_id` between the two sides, or select residues
by either ID from an atom-spec.

## Command

```
swapasym [structures]  [mode  auto|label|auth]
```

| Option | Meaning |
|---|---|
| `structures` | Atom-spec of structures to operate on. Defaults to all open atomic models. |
| `mode auto` | Toggle: auth → label on first run, label → auth on next run (default). Raises when the structure is in a mixed state. |
| `mode label` | Force `chain_id := label_asym_id`. |
| `mode auth` | Force `chain_id := auth_asym_id` (original PDB chain). |

After the first call the bundle attaches two custom Residue attributes
usable from atom-specs:

```
select ::auth_asym_id="A"
select ::label_asym_id="E"
```

## Example

```
open 4hhb
info chains           # 4 polymer chains: A B C D (auth)
swapasym              # 14 unique chain_ids: A-D polymer,
                      #                      E-J HEM/PO4,
                      #                      K-N waters
select /E,F,G,H,I,J   # select all HEM / PO4 via standard atom-spec
swapasym              # back to 4 chains A B C D
```

Structures that were not loaded from mmCIF (plain `.pdb` files) have no
`mmcif_chain_id` information, so `swapasym` raises a UserError when you
try to use it. Reload the structure from a `.cif` file to proceed.

## Installation

### From source (development)

Using [echidna](https://github.com/N283T/echidna) (`echi`), a thin CLI for
ChimeraX bundle development:

```bash
echi build
echi install
```

Or without the helper, using ChimeraX directly:

```bash
ChimeraX --nogui --exit --cmd 'devel install .'
```

### From wheel

```bash
ChimeraX --nogui --exit --cmd 'toolshed install dist/ChimeraX_SwapAsym-0.1.0-py3-none-any.whl'
```

## Development

```bash
echi build                          # Build wheel
echi install                        # Install to ChimeraX
echi run --script scripts/smoke.cxc # Install and run smoke script
```

### Tests

Pytest suite uses stubbed ChimeraX imports and can run under any Python
≥3.11 without a ChimeraX installation:

```bash
uv run --with pytest pytest
```

## License

MIT
