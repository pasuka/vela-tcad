# Reference TCAD CSV Fixtures

This directory is reserved for neutral, text-based reference TCAD exports used
for cross-checking Vela decks. Public paths and scripts intentionally avoid
commercial tool names.

Each device subdirectory can hold:

- `nodes.csv` with `id,x_um,y_um`
- `elements.csv` with `id,node0,node1,node2,region,material`
- `contacts.csv` with `name,node_ids,region`
- `doping.csv` with `node_id,donors_cm3,acceptors_cm3`

Use `scripts/convert_tcad_export.py` to generate Vela `unit_scaling` decks.

Checked-in validation chains currently cover:

- `pn_diode`
- `nmos2d`
- `pmos2d`
- `ldmos2d`
- `igbt2d`
