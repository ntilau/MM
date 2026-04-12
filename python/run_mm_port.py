from __future__ import annotations

import argparse
import sys

from mm_port.scripts import (
    BifurcationE,
    BifurcationH,
    HildebrandFull,
    HildebrandHalf,
    HildebrandSemiAuto,
    Riblet,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run pure-Python MM scripts")
    script_aliases = {
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
    parser.add_argument("script", help="Script to run")
    parser.add_argument("--no-plot", action="store_true", help="Disable plotting")
    args = parser.parse_args()

    script_key = args.script.lower()
    if script_key not in script_aliases:
        allowed = ", ".join(sorted(script_aliases))
        print(f"Unknown script '{args.script}'. Allowed values: {allowed}", file=sys.stderr)
        raise SystemExit(2)

    plot = not args.no_plot
    script = script_aliases[script_key]
    if script == "BifurcationE":
        sf, sinfo, err = BifurcationE(plot=plot)
    elif script == "BifurcationH":
        sf, sinfo, err = BifurcationH(plot=plot)
    elif script == "HildebrandHalf":
        sf, sinfo, err = HildebrandHalf(plot=plot)
    elif script == "HildebrandSemiAuto":
        sf, sinfo, err = HildebrandSemiAuto(plot=plot)
    elif script == "HildebrandFull":
        sf, sinfo, err = HildebrandFull(plot=plot)
    elif script == "Riblet":
        sf, sinfo, err = Riblet(plot=plot)
    else:
        sf, sinfo, err = BifurcationE(plot=plot)

    print(f"Computed {len(sf)} frequency points, {len(sinfo)} ports, fatal={err.get('fatal')}")


if __name__ == "__main__":
    main()
