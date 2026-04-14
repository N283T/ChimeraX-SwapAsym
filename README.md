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

This bundle copies both into custom attributes so you can switch the primary
`chain_id` between the two sides, or select residues by either ID from an
atom-spec.

## Command

```
swapasym [structures]  [mode  auto|label|auth]
```

| Option | Meaning |
|---|---|
| `structures` | Atom-spec of structures to operate on. Defaults to all open atomic models. |
| `mode auto` | Toggle: auth → label on first run, label → auth on next run. (default) |
| `mode label` | Force `chain_id := label_asym_id`. |
| `mode auth`  | Force `chain_id := auth_asym_id` (original PDB chain). |

After the first call the bundle attaches two custom Residue attributes
that are usable from atom-specs:

```
select ::auth_asym_id="A"
select ::label_asym_id="E"
```

## Example

```
open 1fat
info chains           # 4 chains: A B C D (auth)
swapasym              # 20 chains: A..T   (label)
swapasym              # back to 4 chains A B C D
```

## Installation

### From source (development)

```bash
echi install
```

### From wheel

```bash
ChimeraX --nogui --exit --cmd 'toolshed install dist/ChimeraX-SwapAsym-0.1.0-py3-none-any.whl'
```

## Development

```bash
echi build                          # Build wheel
echi install                        # Install to ChimeraX
echi run --script scripts/smoke.cxc # Install and run smoke test
```

## License

MIT
