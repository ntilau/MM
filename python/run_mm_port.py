from __future__ import annotations
"""Small CLI for running the Python modal-matching benchmark scripts."""

import argparse
import sys
from typing import Callable

from mm_port.scripts import (
    BifurcationE,
    BifurcationH,
    HildebrandFull,
    HildebrandHalf,
    HildebrandSemiAuto,
    Riblet,
)

ScriptFn = Callable[[bool], tuple[list, list, dict]]

SCRIPT_ALIASES = {
    "bifurcatione": "BifurcationE",
    "bifurcationh": "BifurcationH",
    "hildebrandhalf": "HildebrandHalf",
    "hilderbrandhalf": "HildebrandHalf",
    "hildebrandsemiauto": "HildebrandSemiAuto",
    "hilderbrandsemiauto": "HildebrandSemiAuto",
    "hildebrandfull": "HildebrandFull",
    "hilderbrandfull": "HildebrandFull",
    "riblet": "Riblet",
}

SCRIPT_RUNNERS: dict[str, ScriptFn] = {
    "BifurcationE": BifurcationE,
    "BifurcationH": BifurcationH,
    "HildebrandHalf": HildebrandHalf,
    "HildebrandSemiAuto": HildebrandSemiAuto,
    "HildebrandFull": HildebrandFull,
    "Riblet": Riblet,
}


def main() -> None:
    """Parse CLI arguments, resolve aliases, run one script, and print summary."""
    parser = argparse.ArgumentParser(description="Run pure-Python MM scripts")
    parser.add_argument("script", help="Script to run")
    parser.add_argument("--no-plot", action="store_true", help="Disable plotting")
    args = parser.parse_args()

    script_key = args.script.lower()
    if script_key not in SCRIPT_ALIASES:
        allowed = ", ".join(sorted(SCRIPT_ALIASES))
        print(f"Unknown script '{args.script}'. Allowed values: {allowed}", file=sys.stderr)
        raise SystemExit(2)

    plot = not args.no_plot
    script_name = SCRIPT_ALIASES[script_key]
    script_runner = SCRIPT_RUNNERS[script_name]
    sf, sinfo, err = script_runner(plot=plot)

    print(f"Computed {len(sf)} frequency points, {len(sinfo)} ports, fatal={err.get('fatal')}")


if __name__ == "__main__":
    main()
