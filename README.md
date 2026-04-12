# MM — Modal-Matching Waveguide Solver

A MATLAB toolbox (with a pure-Python port) for computing generalised scattering
matrices (GSMs) of rectangular waveguide structures via the **modal-matching**
(also known as *mode-matching*) method.

## Overview

The solver assembles multi-port waveguide devices from elementary building
blocks called **TwoPortDevices** (a cascade of rectangular waveguide sections)
and **Nto1Connections** (junctions between N guides and a single wider guide).
For each frequency point the GSM of the full structure is obtained by

1. Computing mode indices and eigenfunctions for every cross-section.
2. Evaluating overlap coupling integrals between modes at each interface.
3. Building per-step scattering matrices and cascading them with phase delays.
4. Condensing the global GSM to expose only the user-defined open ports.

Reference designs implemented in the project scripts include the **Riblet** and
**Hildebrand** direct waveguide couplers and two bifurcation benchmarks
(**BifurcationE** and **BifurcationH**) validated against HFSS full-wave data.

---

## Repository layout

```
MM/
├── matlab/               # MATLAB solver library and benchmark scripts
│   ├── lib/              # Reusable solver functions (add to MATLAB path)
│   └── *.m               # Top-level benchmark scripts
├── python/               # Pure-Python implementation of the same API
│   ├── run.py            # CLI entry point
│   ├── core.py           # Core numerical and solver functions
│   ├── scripts.py        # Top-level benchmark scripts
│   └── lib/              # One module per matlab/lib/*.m function
├── HFSSe.csv             # HFSS reference data for BifurcationE
└── HFSSh.csv             # HFSS reference data for BifurcationH
```

---

## MATLAB

### Requirements

- MATLAB R2014b or later.
- No additional toolboxes are required.

### Getting started

All top-level scripts bootstrap the library path themselves. Either run a
script from within its own directory or from the repository root:

```matlab
% From the repository root
run('matlab/BifurcationE.m')
run('matlab/BifurcationH.m')
run('matlab/Riblet.m')
run('matlab/HildebrandHalf.m')
run('matlab/HildebrandSemiAuto.m')
run('matlab/HildebrandFull.m')
```

To use the library from your own script, add `lib/` to the path first:

```matlab
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'lib'));
```

### Benchmark scripts

| Script | Structure | Ports | Frequency |
|---|---|---|---|
| `BifurcationE.m` | E-plane T-junction | 3 | 10–25 GHz |
| `BifurcationH.m` | H-plane T-junction | 3 | 10–25 GHz |
| `Riblet.m` | Riblet broadwall coupler | 4 | 7–15 GHz |
| `HildebrandHalf.m` | Hildebrand coupler (half, exploits symmetry) | 2 | 13–15 GHz |
| `HildebrandSemiAuto.m` | Hildebrand coupler (half, explicit assembly) | 2 | 13–15 GHz |
| `HildebrandFull.m` | Hildebrand coupler (full 4-port, 6 devices) | 4 | 13–15 GHz |

The bifurcation scripts overlay `HFSSe.csv` / `HFSSh.csv` when the files are
present at the repository root.

### Library reference

#### Solver pipeline

| Function | Signature | Description |
|---|---|---|
| `MultiPortDevice` | `(WGS, Nto1, OpenPorts, ConnectedPorts, FS, flag, Options) → [Sf, Sinfo, …, Error]` | Top-level entry point: validates, solves, and optionally draws |
| `MultiPortDeviceValidate` | `(…) → [WGS, Nto1, …, Symmetry, Topology, Error]` | Validates all inputs, fills defaults, builds topology |
| `MultiPortDeviceSolve` | `(…) → [Sf, Sinfo, Error, …]` | Runs frequency sweep and assembles condensed GSM |
| `MultiPortDeviceTopology` | `(…) → Topology` | Builds internal port connectivity graph |
| `MultiPortDeviceDraw` | `(…)` | Renders 3-D exploded device view |

#### Scattering matrix primitives

| Function | Description |
|---|---|
| `SingleStep` | Scattering matrix at a single waveguide cross-section step |
| `SingleCascade` | Cascade two back-to-back steps with a propagation delay |
| `MultiStep` | Multi-section cascade for an entire `TwoPortDevice` arm |
| `Cascade` | Full `WaveGuideStructure` cascade |
| `CondenseGSM` | GSM condensation via the Selleri method |
| `RenormalizeGSM` | Re-normalise a multi-port GSM to new reference impedances |
| `Renormalize` | Re-normalise a two-port S-matrix |

#### Junction handling

| Function | Description |
|---|---|
| `Nto1Junction` | Builds the scattering matrix of an N-to-1 rectangular junction |
| `Nto1DeviceValidate` | Validates an `Nto1Connection` definition |
| `Nto1DeviceDraw` | 3-D visualisation of an Nto1 junction |

#### Modal analysis

| Function | Description |
|---|---|
| `OrderModes` | Assigns TE/TM mode indices for a cross-section respecting symmetry |
| `EigenModes` | Computes $k_x$, $k_y$, $k_z$ eigenvalues for all retained modes |
| `OneModeEigens` | Eigenvalues for a single $(m,n)$ mode pair |
| `NormCoeff` | Mode normalisation constants $A_h$, $A_e$ |
| `OneModeNormCoeff` | Normalisation constant for a single mode |
| `DelayMatrix` | Diagonal phase-delay matrix $D = \mathrm{diag}(e^{-j k_z l})$ |
| `WaveNumbers` | Propagation constants $k_z$ at a given free-space wavenumber $k_0$ |
| `DgammaMatrices` | Diagonal $\Gamma$ matrices for a cross-section step |
| `UMatrices` | Unit-normalised coupling matrices |
| `MxxMatrices` | Mode-coupling integral matrices $M_{hh}$, $M_{he}$, $M_{eh}$, $M_{ee}$ |
| `Integrals` | Overlap integrals between two rectangular waveguide modes |

#### Port and GSM utilities

| Function | Description |
|---|---|
| `ExtractPortS` | Extract an $(i,j)$ sub-block from a partitioned GSM |
| `InsertPortS` | Insert a sub-block back into a GSM |
| `ExtractSingleS` | Extract one modal $S_{ij}$ coefficient by mode type and index |
| `FrequencySweepValidate` | Validate and complete a `FrequencySweep` struct |

#### Geometry helpers

| Function | Description |
|---|---|
| `TwoPortDeviceValidate` | Validate and fill a `TwoPortDevice` definition |
| `TwoPortDeviceGetPortSegment` | Retrieve the `WaveGuideSegment` at a given port |
| `TwoPortDeviceInsertPortSegment` | Insert or replace a port segment |
| `TwoPortDeviceDraw` | 3-D rendering of a `TwoPortDevice` |
| `ReverseWaveGuideStructure` | Reverse the section order (flip input ↔ output) |
| `WaveGuideSegmentGetBounding` | Axis-aligned bounding box of a segment |
| `WaveGuideSegmentGetCrossSection` | Cross-section parameters |
| `WaveGuideCapDraw` | Draw end-cap patch for visualisation |
| `WaveGuideConnectionCapDraw` | Draw interface cap at a device connection |
| `ShowSegment` | 3-D plot of a single `WaveGuideSegment` |
| `NotInRect` | Geometry test: point outside a set of rectangles |

#### Post-processing and diagnostics

| Function | Description |
|---|---|
| `GSMDraw` | Plot modal S-parameters in dB vs frequency |
| `RelativePhaseDraw` | Plot phase difference between two modal S-parameters |
| `DumpError` | Print an `Error` struct to stdout; returns `halt` flag |

---

### Data structures

#### `WaveGuideSegment`

A struct describing one rectangular waveguide cross-section.

| Field | Type | Description |
|---|---|---|
| `a`, `b` | scalar | Width and height [m] |
| `l` | scalar | Section length [m] |
| `xo`, `yo`, `zo` | scalar | Centre-position offset [m] |
| `Nmodes` | integer | Number of modes to retain |
| `Nh`, `Ne` | integer | Number of TE and TM modes after ordering |
| `mh`, `nh`, `me`, `ne` | vector | Mode index arrays |
| `kh.x`, `kh.y`, `kh.mn` | vector | TE eigenvalues and propagation constants |
| `ke.x`, `ke.y`, `ke.mn` | vector | TM eigenvalues and propagation constants |
| `Ah`, `Ae` | vector | Normalisation coefficients |
| `D` | diagonal matrix | Propagation-delay matrix |

#### `TwoPortDevice`

`WGS{i}` — a cell whose `.D` field is an ordered cell array of
`WaveGuideSegment` structs forming the cascade of sections for device arm `i`.
Port 1 is `WGS{i}.D{1}` and port 2 is `WGS{i}.D{end}`.

#### `Nto1Connection`

Describes the junction between N input `TwoPortDevice` arms and a single wider
guide.

| Field | Description |
|---|---|
| `SideOne{i}.TwoPortDeviceIndex` | Index of the i-th input device |
| `SideOne{i}.TwoPortDevicePort` | Port number (1 or 2) on that device |
| `SideTwo{1}.TwoPortDeviceIndex` | Index of the single output device |
| `SideTwo{1}.TwoPortDevicePort` | Port number on the output device |
| `zo` | z-coordinate of the junction plane [m] |

#### `FrequencySweep`

| Field | Description |
|---|---|
| `start`, `end` | Sweep limits [Hz] |
| `N` | Number of uniformly-spaced points |
| `f` | *(computed)* Frequency vector [Hz] |

#### `Options`

Solver and draw control struct passed to `MultiPortDevice`.

| Field | Values | Description |
|---|---|---|
| `DeviceSymmetry.Use` | `0` / `1` | Enable geometric symmetry to halve solve time |
| `DeviceSymmetry.Side` | `1` / `2` | Which half of the device to solve (port-1 or port-2 side) |
| `Connections` | `0` / `1` | Enable `ConnectedPorts` internal port merging |

#### `Error`

Struct returned by most library functions.  A `fatal` flag indicates an
unrecoverable problem; call `DumpError` to print it and obtain a `halt` flag.

---

### `GSMDraw` — ModeStruct format

`ModeStruct` is a cell array; each entry describes one trace to plot:

```matlab
ModeStruct{k} = {outPort, inPort, outType, outM, outN, inType, inM, inN, label}
```

| Position | Example | Description |
|---|---|---|
| 1 | `1` | Output port index |
| 2 | `2` | Input port index |
| 3 | `'h'` | Output mode family: `'h'` = TE, `'e'` = TM |
| 4 | `1` | Output mode index *m* |
| 5 | `0` | Output mode index *n* |
| 6 | `'h'` | Input mode family |
| 7 | `1` | Input mode index *m* |
| 8 | `0` | Input mode index *n* |
| 9 | `'md'` | Label tag (used in legend) |

`RelativePhaseDraw` takes the same `ModeStruct` format but requires exactly two
entries and plots their phase difference in degrees.

---

### Typical workflow for a new device

```matlab
%% 1. Define waveguide sections (one WGS per arm, one D per section)
WGS{1}.D{1}.a      = 0.01905;  WGS{1}.D{1}.b      = 0.009525;
WGS{1}.D{1}.Nmodes = 12;       WGS{1}.D{1}.l      = 0.01;
WGS{1}.D{1}.xo     = 0;        WGS{1}.D{1}.yo     = 0;
WGS{1}.D{1}.zo     = 0;

WGS{2}.D{1}.a      = 0.01905;  WGS{2}.D{1}.b      = 0.009525;
WGS{2}.D{1}.Nmodes = 12;       WGS{2}.D{1}.l      = 0.01;
WGS{2}.D{1}.xo     = 0;        WGS{2}.D{1}.yo     = 0;
WGS{2}.D{1}.zo     = 0;

%% 2. Define the N-to-1 junction
Nto1{1}.SideOne{1}.TwoPortDeviceIndex = 1;
Nto1{1}.SideOne{1}.TwoPortDevicePort  = 2;
Nto1{1}.SideTwo{1}.TwoPortDeviceIndex = 2;
Nto1{1}.SideTwo{1}.TwoPortDevicePort  = 1;
Nto1{1}.zo = 0;

%% 3. Define open ports and frequency sweep
OpenPorts{1}.TwoPortDeviceIndex = 1;  OpenPorts{1}.TwoPortDevicePort = 1;
OpenPorts{2}.TwoPortDeviceIndex = 2;  OpenPorts{2}.TwoPortDevicePort = 2;
FS.start = 8e9;  FS.end = 12e9;  FS.N = 51;

%% 4. Set options and solve
ConnectedPorts = {};
Options.DeviceSymmetry.Use  = 0;
Options.DeviceSymmetry.Side = 2;
Options.Connections         = 0;

[Sf, Sinfo, WGS, Nto1, ConnectedPorts, FS, Error] = ...
    MultiPortDevice(WGS, Nto1, OpenPorts, ConnectedPorts, FS, 2, Options);

if Error.fatal; return; end

%% 5. Plot — S21 for TE10 input/output
f = FS.f;
ModeStruct = { {2, 1, 'h',1,0, 'h',1,0, 'S21'} };
GSMDraw(f, Sf, Sinfo, ModeStruct, 1);
xlabel('Frequency (Hz)');  ylabel('|S| (dB)');  legend('S21 TE10');
```

---

## Python implementation

A pure-Python implementation of the solver lives in `python/`. It maintains
MATLAB `lib/` API compatibility with matching function names and equivalent
data conventions.

### Requirements

- Python 3.9+
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/) *(optional; required only for plotting)*

No additional installation is required. Run scripts directly from the
repository root.

### Python layout

```
python/
├── run.py                  # CLI entry point
├── core.py                 # Core numerical and solver functions
├── scripts.py              # Top-level benchmark scripts
└── lib/                    # One module per matlab/lib/*.m function
```

### Running the Python CLI

From the repository root:

```bash
python python/run.py <Script> [--no-plot]
```

| Option | Description |
|---|---|
| `<Script>` | One of `BifurcationE`, `BifurcationH`, `Riblet`, `HildebrandHalf`, `HildebrandSemiAuto`, `HildebrandFull` |
| `--no-plot` | Skip plotting for headless or CI runs |

Script names are case-insensitive. Common misspellings such as
`hilderbrandhalf` are accepted as aliases.

Examples:

```bash
# Single run without plot
python python/run.py BifurcationE --no-plot

# Interactive plot
python python/run.py BifurcationH

# Run all benchmarks headlessly
for s in BifurcationE BifurcationH Riblet HildebrandHalf HildebrandSemiAuto HildebrandFull; do
    echo "=== $s ==="
    python python/run.py "$s" --no-plot
done
```

### Python status

| Script | Runs | HFSS overlay |
|---|---|---|
| `BifurcationE` | yes | yes (`HFSSe.csv`) |
| `BifurcationH` | yes | yes (`HFSSh.csv`) |
| `Riblet` | yes | — |
| `HildebrandHalf` | yes | — |
| `HildebrandSemiAuto` | yes | — |
| `HildebrandFull` | yes | — |

Every MATLAB `lib/*.m` function has a Python counterpart in `python/lib/`.
Some advanced drawing helpers in `python/core.py` (notably
`RelativePhaseDraw`) are present but currently raise `NotImplementedError`.

### Python reference overlay behavior

When `HFSSe.csv` / `HFSSh.csv` are present at the repository root, bifurcation
scripts overlay full-wave reference traces on top of MM curves. The overlay
logic in `python/scripts.py`:

- Pairs reference columns to MM traces by RMSE.
- Preserves MM line colors for the paired reference traces.
- Supports optional dB bias correction for `BifurcationH` to improve visual
    parity.
