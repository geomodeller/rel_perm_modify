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
- Post-processing scripts that read grid arrays as `(NX, NY, NZ)` from the SR3
  (e.g., `script/sr3_reader.py` / the RPT-update scripts) only see the
  fundamental grid; refined-block values need extra handling if you want to
  use velocities from inside the LGR.
