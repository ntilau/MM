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
- API coverage:
  - Every MATLAB `lib/*.m` has a Python counterpart module with the same function name.
- Remaining work:
  - Some advanced functions are present but marked `NotImplementedError` in `mm_port/core.py` (mainly drawing/validation and complex multi-topology paths).

## Run

```bash
PYTHONPATH=python python python/run_mm_port.py BifurcationE --no-plot
PYTHONPATH=python python python/run_mm_port.py BifurcationH --no-plot
```
