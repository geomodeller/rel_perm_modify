# How to Apply Local Grid Refinement (LGR) Around a Well in CMG GEM

Applied to: `runs/field_porthos_with-lgr/` (copy of `runs/field_porthos_no-update/`
with 5 m × 5 m horizontal refinement around CO2 injector `5-S1_INJ`).
Deck validated with `gm202420.exe -f datafile.dat -checkonly` → Normal Termination, 0 severe errors.

---

## 1. The `REFINE` keyword (Reservoir Description section)

GEM builds LGR with the `*REFINE` keyword, placed in the **GRID (Reservoir
Description) section** — i.e., after the fundamental grid definition
(`GRID CORNER 38 42 43` + corner arrays) and **before `END-GRID`**:

```
*REFINE  i1(:i2)  j1(:j2)  k1(:k2)  *INTO  ni nj nk
```

Each fundamental (parent) block in the I/J/K range is subdivided into
`ni × nj × nk` child blocks. Children split the parent geometry evenly, so:

- parent areal size 100 m × 100 m (this model, verified from `xcorn.inc`/`ycorn.inc`)
- `INTO 20 20 1` → children of **5 m × 5 m**, layer thickness unchanged (`nk = 1`)

What was added to `datafile.dat` (right after the `zcorn.inc` include):

```
REFINE 19:20 11:12 1:43 INTO 20 20 1
```

The range 19:20 × 11:12 covers every column the well penetrates:
(19,12) K4–15, (20,12) K15–36, (20,11) K36–43.

Other refinement forms in the manual (not used here):

- `*REFINE i j k` on an already-refined block → **nested LGR** (e.g., 100 m → 20 m → 5 m
  in two stages, gentler size contrast at the LGR boundary).
- `*REFINE ... *HYBRID` → hybrid **radial** children inside a Cartesian parent,
  designed for near-well flow (an alternative worth considering for a single injector).

## 2. Property inheritance — and overriding properties inside the LGR

**No property file changes are needed.** Arrays given on the fundamental grid
(`POR`, `PERM*`, `*RTYPE`, corner geometry, etc.) are automatically inherited by
the child blocks from their parent. In this model all children simply inherit
(`porosity.inc`, `permeability.inc`, `rocktype.inc` copied verbatim).

To change properties **inside the refined blocks only**, use the `*RG`
(refined grid) array qualifier:

```
<array keyword> RG <parent address> <array reading option> <values>
```

All of the following were verified to parse cleanly in this deck with
`gm202420.exe -checkonly` (GEM 2024.20):

```
** all 400 children of parent (19,12,5) get porosity 0.30
POR   RG 19 12 5 CON 0.30

** same for permeability
PERMI RG 19 12 5 CON 250.0

** ranges of parent blocks work too
POR   RG 20 12 6:8 CON 0.28

** rock type as well (ROCKFLUID section, after the rocktype.inc include)
RTYPE RG 19 12 5 CON 10
```

Rules and options:

- **Placement / order**: the `RG` line must appear *after* the `REFINE`
  keyword and *after* the base array assignment it overrides (e.g., after the
  `porosity.inc` include for `POR`, after `rocktype.inc` for `RTYPE`).
  Last assignment wins.
- **Per-child values**: replace `CON` with `ALL` and give one value per child
  in local I-fastest order — for `INTO 20 20 1` that is 20×20×1 = 400 values
  per parent, e.g. `POR RG 19 12 5 ALL 400*0.3` or an explicit list. This is
  how each 5 m cell gets its own individual value.
- **What can be overridden**: essentially any block array (`POR`,
  `PERMI/J/K`, `RTYPE`, `NULL`, …). Note `PERMJ EQUALSI` in this deck applies
  globally, so children with an overridden `PERMI` get their J-perm updated
  automatically.
- **RPT-update workflow**: `RTYPE RG ... ALL` is the mechanism needed if
  velocity-based RPT updates should ever act *inside* the LGR — each child can
  carry its own rock type. The catch: `script/sr3_reader.py` only reshapes the
  fundamental 38×42×43 arrays, so reading per-child velocities from the SR3
  needs extra handling first (see §5).

## 3. Well perforations must address child blocks (UBA with `/`)

Once a block is refined, a well can no longer perforate the parent — every
perf UBA must point at a **leaf** (child) block. The UBA syntax appends the
local child address after a slash:

```
parent_i parent_j parent_k / child_i child_j child_k
19 12 4 / 11 11 1
```

Child indices are 1-based, run in the same direction as the fundamental axes
(child 1 at the low-I/low-J face), so with 20 children per 100 m block, child
`n` spans `[5·(n−1), 5·n)` m from the parent's low edge.

### Mapping used for `5-S1_INJ`

The well trajectory (from `LAYERXYZ`) runs through parent-block **centers**:
x = 3850 in block 19 (spans 3800–3900), x = 3950 in block 20, y = 3150 in
j = 12, y = 3050 in j = 11. A block center lies exactly on the face between
children 10 and 11, so the **upper/right child (index 11) was chosen
consistently**; the ≤ 2.5 m offset is negligible at field scale. (If you want
the well exactly at a child center, use an odd subdivision such as
`INTO 21 21 1` → 4.76 m cells.)

| Trajectory section | Parent UBA | Child UBA |
|---|---|---|
| Vertical, K = 4–15 | `19 12 K` | `19 12 K / 11 11 1` |
| Lateral jog at K = 15 (x 3850→3950) | `19 12 15`, `20 12 15` | 20 rows: `19 12 15 / 11..20 11 1`, `20 12 15 / 1..10 11 1` |
| Vertical, K = 16–36 | `20 12 K` | `20 12 K / 11 11 1` |
| Lateral jog at K = 36 (y 3150→3050) | `20 12 36`, `20 11 36` | 20 rows: `20 12 36 / 11 10..1 1`, `20 11 36 / 11 20..11 1` |
| Vertical, K = 37–43 | `20 11 K` | `20 11 K / 11 11 1` |

The two 50 m lateral segments were split into 10 sub-segments of 5 m each
(one per traversed child), with the `LAYERXYZ` entry/exit coordinates linearly
interpolated and the segment length divided proportionally. The `FLOW-FROM`
chain was renumbered to stay sequential. Result: 39 original perf rows → **75
child-block perf rows**, total perforated length preserved exactly
(269.8712 m), every child UBA unique.

## 4. What changed in `runs/field_porthos_with-lgr/`

| File | Change |
|---|---|
| `datafile.dat` | `REFINE 19:20 11:12 1:43 INTO 20 20 1` inserted in GRID section |
| `well_config.inc` | `PERF GEOA` + `LAYERXYZ` tables rewritten to 75 child-block UBAs |
| all other `.inc` | copied unchanged (properties inherited by children) |

## 5. Model size and practical caveats

- Block count: 68,628 fundamental + 172 × 400 children = **137,428 blocks**
  (GEM check run confirms). Expect a substantially longer runtime, and smaller
  time steps early on (5 m cells next to an injector throttle the CFL/throughput
  limit).
- GEM warns that some layer-1 parents in the refined range are null / have no
  thickness (e.g., block (19,11,1)) — benign; null parents are skipped.
- The 20:1 size jump at the LGR boundary is legal but numerically abrupt; if
  convergence suffers, use nested refinement (100 → 20 → 5 m) or shrink `INTO`.
- ~~Post-processing scripts that read grid arrays as `(NX, NY, NZ)` from the SR3
  only see the fundamental grid.~~ Resolved — see §6: `sr3_reader.py` is now
  LGR-aware.

## 6. Reading refined-grid results from the SR3 (verified layout)

A short LGR test run (GEM 2024.20) was used to decode how the SR3 stores
refined grids. Layout, all verified empirically:

- **Complete storage** = fundamental cells first (68,628, natural I-fastest
  order — refined parents keep their slot but become inactive), then all child
  cells appended (172 grids × 400 = 68,800 → 137,428 total).
- **Refined grids are created I-fastest, then J, then K** over the `REFINE`
  range: (19,11,1), (20,11,1), (19,12,1), (20,12,1), (19,11,2), … , (20,12,43).
- **Within a refined grid, children are stored local-I-fastest** (confirmed by
  matching child depth gradients against the fundamental structural dip) —
  the same ordering `*RG ... *ALL` expects on input.
- The mapping lives in the SR3 `GRID` datasets: `ICSTCG` (per cell, the number
  of the refined grid it hosts; 0 = none), `IGNTID/JD/KD` (per-grid dims),
  `IGNTNC` (cumulative cell counts → each grid's slice of complete storage),
  `ICSTPS` (complete → packed/active mapping, used to unpack property arrays).
- Every grid property in `OUTSRF GRID` (PRES, SG, VELOCRC, …) is written for
  child cells too; per-perforation well rates appear in
  `TimeSeries/LAYERS` (75 layers = the 75 child-block perfs).

### sr3_reader.py support

- `get_grid()` was fixed to slice geometry arrays to the fundamental grid
  (it crashed on LGR SR3s before) and now sets `sr3.grid.n_grids` /
  `sr3.grid.n_fund_cells`.
- **`get_lgr_table(sr3)`** → one record per refined grid:
  `{'grid', 'parent_ijk', 'dims', 'start', 'stop'}` where start/stop slice the
  child section of complete storage.
- **`get_grid_properties_lgr(sr3, NX, NY, NZ, nt)`** → `(fund, children,
  lgr_table)`: `fund[name]` is the usual `(nt, NZ, NY, NX)` array;
  `children[name]` is `(nt, n_child_cells)` — slice with
  `lgr_table[m]['start']:['stop']` and reshape `(20, 20)` (local j, local i)
  for one parent. Backward compatible: on a non-LGR SR3, `children` is empty
  and `fund` matches `get_grid_properties()` exactly.

Example — CO2 gas velocity magnitude in the 5 m cells at the last time step:

```python
sr3 = read_SR3('datafile.sr3')
nt  = len(sr3.times['Days'])
fund, ch, tab = get_grid_properties_lgr(sr3, 38, 42, 43, nt)
vel = np.sqrt(ch['VELGXRC'][-1]**2 + ch['VELGYRC'][-1]**2 + ch['VELGZRC'][-1]**2)
rec = next(r for r in tab if r['parent_ijk'] == (19, 12, 8))
vel_parent = vel[rec['start']:rec['stop']].reshape(20, 20)   # [local j, local i]
```

(Sanity check from the test run: max child velocity 24.6 m/day lands in parent
(19,12,8) child (10,11) — directly beside the well child (11,11).)

## 7. RPT-update workflow with LGR

- **`runs/field_porthos_with-update_lgr/`** = `runs/field_porthos_with-update/`
  + the `REFINE` line in `template.dat` + the child-UBA `well_config.inc`.
- **`script/run_cmg_with_rpt_update_4_field_lgr.py`** runs it: same monthly
  restart loop as the field script (30 RPT bins of 1.016 m/day, `**$`
  placeholder activation), but velocities are read for fundamental + child
  cells, and `rocktype.inc` is written as `*RTYPE *ALL` (fundamental) followed
  by one `*RTYPE *RG i j k *ALL` block per refined parent (400 values,
  local-I-fastest — the SR3 child order can be written back verbatim).
  A full RG-style `rocktype.inc` (172 blocks) passed `-checkonly` cleanly.
