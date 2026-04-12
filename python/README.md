# Pure Python Port (MM)

This folder contains a pure-Python port of the modal-matching solver API for the repository.

## Layout

- `mm_port/core.py`: core numerical and solver functions.
- `mm_port/lib/*.py`: one Python module per MATLAB file in `lib/`, preserving function names.
- `mm_port/scripts.py`: top-level script equivalents.
- `run_mm_port.py`: CLI for running ported scripts.

## Status

- Fully runnable in pure Python:
  - `BifurcationE`
  - `BifurcationH`
  - `Riblet`
  - `HildebrandHalf`
  - `HildebrandSemiAuto`
  - `HildebrandFull`
- Plot parity with MATLAB benchmarks:
  - Both bifurcation scripts can overlay `HFSSe.csv` / `HFSSh.csv` when plotting.
- API coverage:
  - Every MATLAB `lib/*.m` has a Python counterpart module with the same function name.
- Remaining work:
  - Some advanced functions are present but marked `NotImplementedError` in `mm_port/core.py` (notably `RelativePhaseDraw`).

## Run

```bash
python python/run_mm_port.py BifurcationE --no-plot
python python/run_mm_port.py BifurcationH --no-plot
python python/run_mm_port.py Riblet --no-plot
python python/run_mm_port.py HildebrandHalf --no-plot
python python/run_mm_port.py HildebrandSemiAuto --no-plot
python python/run_mm_port.py HildebrandFull --no-plot
```

Run all structures in one pass:

```bash
for s in BifurcationE BifurcationH Riblet HildebrandHalf HildebrandSemiAuto HildebrandFull; do
  echo "=== $s ==="
  python python/run_mm_port.py "$s" --no-plot
done
```

Launch an interactive plot window:

```bash
python python/run_mm_port.py BifurcationH
```

## Reference Overlay Notes

- Bifurcation plots can overlay HFSS reference CSV traces.
- The overlay logic in `mm_port/scripts.py` preserves MM line colors for references.
- For `BifurcationH`, reference curves are globally paired to MM traces and can use affine dB alignment for better visual parity.
