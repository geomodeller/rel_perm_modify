# How to Compute Relative Permeability in CMG GEM 2024.20

> Source: CMG GEM 2024.20 manual (`./GEM/Content/GEM/`)
> Relevant pages: `Rock Fluid Properties/SWT_*`, `Rock Fluid Properties/SLT_*`, `Rock Fluid Properties/KROIL_*`, `Rock Fluid Properties/KRGAS_*`, `Rock Fluid Properties/NotesOnRockFluidProperties.htm`, `Fluid Model/PHASEID_*`, `Fluid Model/FLASH-METHOD_*`, `Fluid Model/NotesOnFluidModel.htm`

---

## 1. Key Concept: Phase-Based, Not Component-Based

GEM assigns relative permeability to **phases** (gas, oil, water), not to individual components.
Each component's flux = (phase Darcy flux) × (mole fraction of that component in the phase).

GEM supports at most **two hydrocarbon phases** (one gas + one oil) simultaneously, plus one aqueous phase.

---

## 2. Input Keywords for Relative Permeability Tables

### `*SWT` — Water-Oil Relative Permeability Table

**Water-wet columns:** `Sw | krw | krow | (Pcow) | (Pcowi)`  
**Oil-wet columns:** `So | kro | krwo | (Pcwo) | (Pcwoi)`

- Must appear after `*RPT` in the Rock Fluid Properties section
- One `*SWT` required per rock type (`*RPT`)
- Entries in **increasing Sw (or So)** order; first entry = connate saturation
- Capillary pressure columns optional (default = 0); both drainage and imbibition needed for hysteresis
- `Pcow = Po − Pw` (water-wet); `Pcwo = Pw − Po` (oil-wet)

### `*SLT` / `*SGT` — Liquid-Gas Relative Permeability Table

**Water-wet columns:** `Sl (or Sg) | krg | krog | (Pcog) | (Pcogi)`  
**Oil-wet columns:** `Sl (or Sg) | krg | krwg | (Pcog) | (Pcogi)`

- `*SLT`: indexed by liquid saturation `Sl = So + Swcon`; entries in increasing Sl order; last entry Sl = 1
- `*SGT`: indexed by gas saturation Sg; entries in increasing Sg order; first entry Sg = 0
- Only one of `*SLT` or `*SGT` per `*RPT`
- **Consistency requirement:** `krog` at Sl=1 (or Sg=0) must equal `krow` at connate water in `*SWT` for the same rock type
- `Pcog = Pg − Po`; decreases monotonically with increasing Sl

### Column Definitions

| Symbol | Meaning |
|--------|---------|
| `krw` | kr to water at given Sw |
| `krow` | kr to oil in presence of water (2-phase oil-water, no gas) — function of Sw |
| `krg` | kr to gas at given Sl or Sg |
| `krog` | kr to oil in presence of gas and connate water (2-phase gas-oil) — function of Sg |
| `Pcow` | Oil-water capillary pressure (drainage) |
| `Pcowi` | Oil-water capillary pressure (imbibition); if omitted = Pcow, no hysteresis |

### `*SMOOTHEND` — Endpoint Smoothing (applies to both `*SWT` and `*SLT`)

Controls interpolation where kr transitions from zero to non-zero:

| Option | Behavior |
|--------|---------|
| `*QUADGEM` | **Default in GEM**; recommended |
| `*QUAD` | Quadratic, similar to IMEX/STARS |
| `*CUBIC` | Cubic interpolation |
| `*LINEAR` | Linear (was default before GEM 2014.11) |
| `*POWERQ m e` | Power-law; falls back to `*QUAD` if exponent outside [1.5, 4.0] |
| `*POWERC m e` | Power-law; falls back to `*CUBIC` if outside [1.5, 4.0] |
| `*ON` | Equivalent to `*QUADGEM` |
| `*OFF` | Equivalent to `*LINEAR` |

---

## 3. Three-Phase kro — `*KROIL`

When all three phases (gas, oil, water) are present, GEM cannot use a single two-phase table for oil. It combines `krow` and `krog` via a three-phase model.

**Default:** `*KROIL *STONE2 *SWSG`

### Stone's 2nd Model (Aziz & Settari normalization)

```
kro = krocw × ((krow/krocw + krw) × (krog/krocw + krg) − krw − krg)
```

where:
- `krocw` = oil kr at connate water (first `krow` value in `*SWT`)
- `krow(Sw)` — looked up from `*SWT`
- `krog(Sg)` — looked up from `*SLT`
- `krw(Sw)` — looked up from `*SWT`
- `krg(Sg)` — looked up from `*SLT`
- Result clipped to 0 if negative

### Other Available Models

| Keyword | Model |
|---------|-------|
| `*KROIL *STONE2 *SWSG` | Stone's 2nd, krow(Sw) and krog(Sg) — **default** |
| `*KROIL *STONE2 *SO` | Stone's 2nd, krow and krog as functions of So |
| `*KROIL *STONE1 *SWSG` | Stone's 1st (Fayers & Matthews Som) |
| `*KROIL *LINEAR_ISOPERM` | Linear Isoperm method (Baker 1988) |
| `*KROIL *SEGREGATED` | Segregated model (Baker 1988) |

### `*KRGAS` — Single Hydrocarbon Phase

When only one hydrocarbon phase is present and identified as gas:

| Keyword | Behavior |
|---------|---------|
| `*KRGAS *KRG` | Use `krg(Sg)` from `*SLT` — **default** |
| `*KRGAS *KRO` | Use `krow(Sw)` from `*SWT` (legacy behavior) |

---

## 4. Application to CO2 + CH4 + Brine System

### 4.1 Phase Determination — Automatic via EOS

GEM determines 2-phase vs 3-phase **automatically** at every gridblock every timestep. The user does not specify this.

#### Hydrocarbon Phase Count (EOS Stability Test)

At each Newton iteration, GEM checks the Tangent Plane Distance (TPD):

- **TPD ≥ 0** → single hydrocarbon phase (CO2 + CH4 fully miscible) → **2-phase total** (hydrocarbon + brine)
- **TPD < 0** → mixture unstable → flash calculates two hydrocarbon phases (CO2-rich dense "oil" + CH4-rich gas) → **3-phase total**

The phase boundary is controlled by EOS parameters:
- `*EOSTYPE` — PR or SRK
- `*PCRIT`, `*TCRIT`, `*VCRIT` — critical properties
- `*BIN` — binary interaction coefficients (most sensitive for CO2-CH4)

#### Aqueous Phase — Always Separate

Brine is always tracked independently from the hydrocarbon EOS.

#### Single-Phase Identity — `*PHASEID`

When the EOS finds only one hydrocarbon phase, GEM must label it "gas" or "oil":

| Keyword | Method |
|---------|--------|
| `*PHASEID *DEN` | Compare density to reference gas/oil densities — **default** |
| `*PHASEID *CRIT` | Compare molar volume to pseudo-critical volume (Gosset et al.) |
| `*PHASEID *TCMIX` | Compare block T to pseudo-critical T (Li mixing rule) |
| `*PHASEID *GAS` | Force all single-phase blocks to gas |
| `*PHASEID *OIL` | Force all single-phase blocks to oil |
| `*REFDEN dvref` | Compare density to user-specified reference |

#### CO2 Dissolution in Brine — OGW Flash

When using `*SOLUBILITY` or `*HENRY-CORR-XXX` keywords (Henry's law for CO2 dissolution):
- `*FLASH-METHOD-OGW` is **automatically activated**
- Enables special Oil/Gas/Water three-phase flash coupling hydrocarbon EOS with aqueous solubility

### 4.2 kr Assignment by Phase

#### Case 1 — 2-Phase (CO2+CH4 gas + brine)

| Phase | Components | kr source |
|-------|-----------|-----------|
| Gas | CO2 + CH4 (shared) | `krg` from `*SLT` at Sl = Sw |
| Aqueous | Brine | `krw` from `*SWT` at Sw |

CO2 and CH4 share the same `krg`. Individual component fluxes are weighted by mole fraction in the gas phase.

#### Case 2 — 3-Phase (CO2-dense "oil" + CH4 gas + brine)

| Phase | Components | kr source |
|-------|-----------|-----------|
| Gas | CH4-rich | `krg` from `*SLT` at Sl = So + Swcon |
| Oil (dense) | CO2-rich | `kro` via Stone's 2nd model |
| Aqueous | Brine | `krw` from `*SWT` at Sw |

`kro` (CO2-rich phase) computed as:
```
kro = krocw × ((krow(Sw)/krocw + krw(Sw)) × (krog(Sg)/krocw + krg(Sg)) − krw(Sw) − krg(Sg))
```

---

## 5. Worked Numerical Example

### Table Setup

```
*ROCKFLUID
*RPT 1

*SWT
**  Sw     krw    krow
    0.20   0.000  1.000    ! connate water Swc=0.20, krocw=1.0
    0.40   0.033  0.667
    0.60   0.133  0.333
    0.80   0.300  0.000
    1.00   1.000  0.000

*SLT
**  Sl     krg    krog
    0.20   0.984  0.000    ! Sl_min = Swc
    0.40   0.870  0.000
    0.60   0.410  0.021
    0.70   0.190  0.090
    0.80   0.075  0.350
    0.90   0.000  0.997
    1.00   0.000  1.000    ! krog at Sl=1 = krow at Swc = 1.0 (must match)

*KROIL *STONE2 *SWSG      ! default
*KRGAS *KRG               ! default
```

### Case 1 — 2-Phase: Sw=0.50, Sg=0.50

- `Sl = Sw = 0.50` (no oil phase)
- `krg(Sl=0.50)` → interpolate SLT between Sl=0.40 (0.870) and Sl=0.60 (0.410): **krg = 0.640**
- `krw(Sw=0.50)` → interpolate SWT between Sw=0.40 (0.033) and Sw=0.60 (0.133): **krw = 0.083**

| Phase | kr |
|-------|----|
| Gas (CO2+CH4) | **0.640** |
| Brine | **0.083** |

### Case 2 — 3-Phase: Sw=0.40, So=0.50 (CO2), Sg=0.10 (CH4)

- `krw(Sw=0.40)` = 0.033; `krow(Sw=0.40)` = 0.667
- `Sl = So + Swc = 0.50 + 0.20 = 0.70`
- `krg(Sl=0.70)` = 0.190 (Sg=0.10 → below critical gas saturation, effectively 0)
- `krog(Sl=0.90, Sg=0.10)` = 0.997

Stone's 2nd:
```
kro = 1.0 × ((0.667/1.0 + 0.033) × (0.997/1.0 + 0.000) − 0.033 − 0.000)
    = (0.700 × 0.997) − 0.033
    = 0.698 − 0.033
    = 0.665
```

| Phase | kr |
|-------|----|
| Gas (CH4) | **0.000** (below critical Sg) |
| Dense CO2 "oil" | **0.665** |
| Brine | **0.033** |

---

## 6. Decision Flowchart

```
Each gridblock, each timestep:
│
├─ Aqueous phase: always present (brine) → krw from *SWT at Sw
│
└─ EOS stability test on hydrocarbon mixture
       │
       ├─ Stable (TPD ≥ 0) → 1 hydrocarbon phase
       │       │
       │       └─ *PHASEID labels it gas or oil
       │           ├─ Gas → krg from *SLT at Sl=Sw   [2-phase total]
       │           └─ Oil → krow from *SWT at Sw      [2-phase total]
       │
       └─ Unstable (TPD < 0) → flash → 2 hydrocarbon phases
               ├─ Gas phase  → krg from *SLT at Sl = So + Swcon  [3-phase total]
               ├─ Oil phase  → kro via *KROIL (Stone's 2nd)
               │               uses krow(Sw) + krog(Sg) + krw(Sw) + krg(Sg)
               └─ Aqueous    → krw from *SWT at Sw
```

---

## 7. Important Notes

1. **Consistency between SWT and SLT:** `krog` at Sl=1 (SLT last row) must equal `krow` at connate water (SWT first row) for the same rock type.
2. **Stone's model can produce negative kro** — GEM clips to 0. This means the oil phase has zero mobility at that saturation state.
3. **`*SMOOTHEND *QUADGEM` is the default** for both `*SWT` and `*SLT` since GEM 2014.11. Use `*LINEAR` for backward compatibility with pre-2014 datasets.
4. **Fluid Model units are independent of `*INUNIT`** — always check the Fluid Model section for correct units (e.g., critical pressure in atm regardless of SI/FIELD setting).
5. **Oil-wet systems:** Do not use `*KROIL` sub-keywords `*SO`/`*SWSG` or `*KRGAS` for oil-wet reservoirs. For oil-wet, `*SWT` first column is So (not Sw) and `*SLT` third column is `krwg` (not `krog`).
6. **`*FLASH-METHOD-OGW`** is automatically enabled when `*SOLUBILITY` or `*HENRY-CORR-XXX` is used — required for proper CO2 dissolution in brine modeling.
