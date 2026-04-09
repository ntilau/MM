# MM — Modal-Matching Waveguide Solver

A MATLAB toolbox for computing generalised scattering matrices (GSMs) of
rectangular waveguide structures via the **modal-matching** (also known as
*mode-matching*) method.

## Overview

The solver assembles multi-port waveguide devices from elementary building
blocks called **TwoPortDevices** (a series of rectangular waveguide sections)
and **Nto1Connections** (junctions between N guides and a single wider guide).
For each frequency point the GSM of the full structure is obtained by

1. Computing mode indices and eigenfunctions for every cross-section.
2. Evaluating coupling integrals between modes at each step.
3. Building per-step scattering matrices and cascading them.
4. Condensing the global GSM to expose only the user-defined open ports.

Reference designs implemented in the project scripts include the
**Riblet** and **Hildebrand** direct waveguide couplers and the two
bifurcation benchmark scripts **BifurcationE** and **BifurcationH**.

---

## Repository layout

```
MM/
├── lib/                  # Library of reusable MATLAB functions
│   ├── MultiPortDevice.m         # Top-level solver entry point
│   ├── MultiPortDeviceSolve.m    # Frequency-sweep engine
│   ├── MultiPortDeviceValidate.m # Input validation
│   ├── MultiPortDeviceTopology.m # Port topology builder
│   ├── MultiPortDeviceDraw.m     # 3-D device visualisation
│   ├── Nto1Junction.m            # N-to-1 junction scattering matrix
│   ├── Nto1DeviceValidate.m      # Nto1 validation
│   ├── Nto1DeviceDraw.m          # Nto1 visualisation
│   ├── MultiStep.m               # Multi-section cascaded scattering matrix
│   ├── SingleStep.m              # Single step scattering matrix
│   ├── SingleCascade.m           # Two-step cascade with delay
│   ├── Cascade.m                 # Full WaveGuideStructure cascade
│   ├── CondenseGSM.m             # GSM condensation (Selleri method)
│   ├── RenormalizeGSM.m          # GSM re-normalisation
│   ├── Renormalize.m             # S-matrix re-normalisation
│   ├── OrderModes.m              # Assign TE/TM mode indices per symmetry
│   ├── EigenModes.m              # Compute kx, ky, kz eigenvalues
│   ├── DelayMatrix.m             # Build diagonal phase-delay matrix
│   ├── NormCoeff.m               # Compute mode normalisation constants
│   ├── OneModeEigens.m           # Eigenvalues for a single (m,n) mode
│   ├── OneModeNormCoeff.m        # Norm coefficient for a single mode
│   ├── DgammaMatrices.m          # Diagonal γ matrices for a step
│   ├── UMatrices.m               # U (unit-normalised) matrices for a step
│   ├── MxxMatrices.m             # Mode-coupling integral matrices (Mhh/Mhe/Meh/Mee)
│   ├── Integrals.m               # Overlap integrals between two modes
│   ├── ExtractPortS.m            # Extract sub-block from a GSM
│   ├── InsertPortS.m             # Insert sub-block into a GSM
│   ├── ExtractSingleS.m          # Extract a single S-parameter by mode
│   ├── FrequencySweepValidate.m  # Newer frequency sweep validator
│   ├── DumpError.m               # Print error structure to stdout
│   ├── GSMDraw.m                 # Plot GSM S-parameters vs frequency
│   ├── RelativePhaseDraw.m       # Plot relative phase between two ports
│   ├── ShowSegment.m             # 3-D plot of a single waveguide segment
│   ├── WaveGuideCapDraw.m        # Low-level cap patch drawing
│   ├── WaveGuideConnectionCapDraw.m # Cap at a connection interface
│   ├── WaveGuideSegmentGetBounding.m    # Bounding-box of a segment
│   ├── WaveGuideSegmentGetCrossSection.m # Cross-section parameters
│   ├── TwoPortDeviceDraw.m       # 3-D plot of a TwoPortDevice
│   ├── TwoPortDeviceGetPortSegment.m    # Retrieve port segment from device
│   ├── TwoPortDeviceInsertPortSegment.m # Insert/replace a port segment
│   ├── TwoPortDeviceValidate.m   # Validate and fill TwoPortDevice
│   ├── ReverseWaveGuideStructure.m # Reverse segment order in a structure
│   ├── NotInRect.m               # Geometry: point outside rectangles test
│   └── WaveNumbers.m             # Compute propagation constants at k0
├── BifurcationE.m        # E-plane bifurcation benchmark and plot script
├── BifurcationH.m        # H-plane bifurcation benchmark and plot script
├── HildebrandFull.m      # Full Hildebrand 4-port coupler
├── HildebrandHalf.m      # Half-coupler (symmetry exploit)
├── HildebrandSemiAuto.m  # Semi-automated Hildebrand assembly
├── Riblet.m              # Riblet direct coupler
├── BifurcationE.mat      # Reference full-wave data for BifurcationE
├── BifurcationH.mat      # Reference full-wave data for BifurcationH
└── old/                  # Archived older versions (not maintained)
```

---

## Data structures

### `WaveGuideSegment`
A struct describing a single rectangular waveguide cross-section and its
electromagnetic modes.

| Field | Description |
|---|---|
| `a`, `b` | Width and height [m] |
| `l` | Section length [m] |
| `xo`, `yo`, `zo` | Centre-position offset [m] |
| `Nmodes` | Number of modes to retain |
| `Nh`, `Ne` | Number of TE and TM modes after ordering |
| `mh`, `nh`, `me`, `ne` | Mode index arrays |
| `kh.x`, `kh.y`, `kh.mn` | TE eigenvalues and propagation constants |
| `ke.x`, `ke.y`, `ke.mn` | TM eigenvalues and propagation constants |
| `Ah`, `Ae` | Normalisation coefficients |
| `D` | Diagonal delay matrix |

### `TwoPortDevice`
A cell entry `WGS{i}` whose `.D` field is a cell array of
`WaveGuideSegment` structs forming the cascade of sections for that arm.

### `Nto1Connection`
Describes the junction between multiple TwoPortDevices.

| Field | Description |
|---|---|
| `SideOne{i}.TwoPortDeviceIndex` | Index of the i-th input device |
| `SideOne{i}.TwoPortDevicePort` | Port (1 or 2) on that device |
| `SideTwo{i}.TwoPortDeviceIndex` | Index of the output device |
| `SideTwo{i}.TwoPortDevicePort` | Port on the output device |
| `zo` | z-position of the junction plane [m] |

### `FrequencySweep`
| Field | Description |
|---|---|
| `start`, `end` | Sweep limits [Hz] |
| `N` | Number of uniformly-spaced points |
| `f` | (computed) Frequency vector [Hz] |

### `Symmetry`
| Field | Description |
|---|---|
| `x`, `y` | 1 = symmetric along x / y wall |
| `H`, `E` | 1 = magnetic / electric wall symmetry |

### `Error`
Struct returned by most functions.  A `fatal` field indicates an
unrecoverable problem; `DumpError` prints the content and returns a
`halt` flag.

---

## Getting started

1. Open MATLAB and set the working directory to the repository root, or
   run any project script directly — the first lines automatically add
   `lib/` to the MATLAB path:
   ```matlab
   projectRoot = fileparts(mfilename('fullpath'));
   addpath(fullfile(projectRoot, 'lib'));
   ```
2. Run any of the top-level scripts, e.g.:
   ```matlab
   run HildebrandFull   % Hildebrand coupler, full structure
   run Riblet           % Riblet coupler
   run BifurcationE    % E-plane bifurcation benchmark
   run BifurcationH    % H-plane bifurcation benchmark
   ```
3. Results are plotted automatically. The bifurcation scripts optionally
   overlay benchmark data loaded from `BifurcationE.mat` and
   `BifurcationH.mat`; the loaded variables remain `HFSS` and `HFSSh`
   inside MATLAB, matching the historical plotting code.

---

## Typical workflow for a new device

```matlab
%% 1. Define waveguide sections
WGS{1}.D{1}.a = 0.01905;  WGS{1}.D{1}.b = 0.009525;
WGS{1}.D{1}.Nmodes = 12;  WGS{1}.D{1}.l = 0.01;
WGS{1}.D{1}.xo = 0;       WGS{1}.D{1}.yo = 0;  WGS{1}.D{1}.zo = 0;
% ... repeat for each device ...

%% 2. Define junctions
Nto1{1}.SideOne{1}.TwoPortDeviceIndex = 1;
Nto1{1}.SideOne{1}.TwoPortDevicePort  = 2;
Nto1{1}.SideTwo{1}.TwoPortDeviceIndex = 2;
Nto1{1}.SideTwo{1}.TwoPortDevicePort  = 1;
Nto1{1}.zo = 0;

%% 3. Define open ports and frequency sweep
OpenPorts{1}.TwoPortDeviceIndex = 1;  OpenPorts{1}.TwoPortDevicePort = 1;
OpenPorts{2}.TwoPortDeviceIndex = 2;  OpenPorts{2}.TwoPortDevicePort = 2;
FS.start = 8e9;  FS.end = 12e9;  FS.N = 51;

%% 4. Solve
ConnectedPorts = {};
Options.DeviceSymmetry.Use = 0;
Options.Connections = 0;
[Sf, Sinfo, WGS, Nto1, ConnectedPorts, FS, Error] = ...
    MultiPortDevice(WGS, Nto1, OpenPorts, ConnectedPorts, FS, 2, Options);

%% 5. Plot
f = FS.f;
ModeStruct = {{1,2,'h',1,0,'h',1,0,'md'}};
GSMDraw(f, Sf, Sinfo, ModeStruct, 1);
```

---

## Dependencies

- MATLAB R2014b or later (uses `isfloat`, `1i` imaginary literal).
- No additional toolboxes are required.
