from __future__ import annotations

import argparse

from mm_port.scripts import BifurcationE, BifurcationH, HildebrandHalf, Riblet


def main() -> None:
    parser = argparse.ArgumentParser(description="Run pure-Python MM scripts")
    parser.add_argument(
        "script",
        choices=["BifurcationE", "BifurcationH", "HildebrandHalf", "Riblet"],
        help="Script to run",
    )
    parser.add_argument("--no-plot", action="store_true", help="Disable plotting")
    args = parser.parse_args()

    plot = not args.no_plot
    if args.script == "BifurcationE":
        sf, sinfo, err = BifurcationE(plot=plot)
    elif args.script == "BifurcationH":
        sf, sinfo, err = BifurcationH(plot=plot)
    elif args.script == "HildebrandHalf":
        sf, sinfo, err = HildebrandHalf(plot=plot)
    elif args.script == "Riblet":
        sf, sinfo, err = Riblet(plot=plot)
    else:
        sf, sinfo, err = BifurcationE(plot=plot)

    print(f"Computed {len(sf)} frequency points, {len(sinfo)} ports, fatal={err.get('fatal')}")


if __name__ == "__main__":
    main()
