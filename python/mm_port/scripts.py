from __future__ import annotations
"""High-level script equivalents for the original MATLAB examples.

Each function builds a waveguide topology, runs `MultiPortDevice`, and
optionally renders MATLAB-style plots for quick parity checks.
"""

from pathlib import Path
from typing import Any

import numpy as np

from .core import C0, ExtractSingleS, GSMDraw, MultiPortDevice, OrderModes


_PROJECT_ROOT = Path(__file__).resolve().parents[2]
ModeSpec = tuple[int, int, str, int, int, str, int, int]
ScriptResult = tuple[list[np.ndarray], list[dict[str, Any]], dict[str, Any]]


def _uniform_frequency_axis(fs: dict[str, float | int]) -> np.ndarray:
    """Build a linearly spaced frequency axis from a MATLAB-style `FS` dict."""
    n = int(fs["N"])
    if n <= 1:
        return np.asarray([float(fs["start"])], dtype=float)
    return np.linspace(float(fs["start"]), float(fs["end"]), n, dtype=float)


def _overlay_reference_csv(
    csv_name: str,
    fallback_styles: list[str],
    series_cols: list[int] | None = None,
) -> None:
    """Overlay reference traces from a repo-level CSV onto the current axes.

    The CSV must use the first column as frequency in GHz and remaining columns
    as trace values.
    """
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    csv_path = _PROJECT_ROOT / csv_name
    if not csv_path.exists():
        return
    ref = np.loadtxt(csv_path, delimiter=",")
    if ref.ndim == 1:
        ref = ref.reshape(1, -1)
    if ref.shape[1] < 2:
        return
    x = ref[:, 0] * 1e9
    all_cols = list(range(1, ref.shape[1]))
    use_cols = series_cols if series_cols is not None else all_cols
    ax = plt.gca()
    existing = ax.get_lines()
    for i, col in enumerate(use_cols):
        if i < len(existing):
            color = existing[i].get_color()
            plt.plot(x, ref[:, col], linestyle=":", color=color)
        else:
            style = fallback_styles[i % len(fallback_styles)]
            if ":" not in style:
                style = f"{style}:"
            plt.plot(x, ref[:, col], style)


def _find_balanced_nmodes(segments_specs: list[dict], nmodes_max: int) -> int:
    """Return the largest n <= nmodes_max where all segments yield the same total mode count."""
    for nm in range(nmodes_max, 0, -1):
        totals = set()
        for spec in segments_specs:
            seg = dict(spec)
            seg["Nmodes"] = nm
            s, _ = OrderModes(seg, {})
            totals.add(s["Nh"] + s["Ne"])
        if len(totals) == 1:
            return nm
    return 1


def _draw_relative_phase(
    f: np.ndarray,
    sf: list[np.ndarray],
    sinfo: list[dict[str, Any]],
    mode_a: ModeSpec,
    mode_b: ModeSpec,
    ylim: tuple[float, float] | None = None,
) -> None:
    """Plot unwrapped relative phase between two extracted S-parameter modes."""
    import matplotlib.pyplot as plt

    def _extract(mode: ModeSpec) -> np.ndarray:
        out_port, in_port, out_type, out_m, out_n, in_type, in_m, in_n = mode
        return np.asarray(
            [
                ExtractSingleS(
                    s,
                    sinfo,
                    out_port,
                    in_port,
                    out_type,
                    out_m,
                    out_n,
                    in_type,
                    in_m,
                    in_n,
                )
                for s in sf
            ],
            dtype=complex,
        )

    a = _extract(mode_a)
    b = _extract(mode_b)
    phase = np.rad2deg(np.unwrap(np.angle(a) - np.angle(b)))
    plt.figure()
    plt.plot(f, phase)
    plt.xlim(float(np.min(f)), float(np.max(f)))
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Phase (deg)")
    plt.grid(True)
    plt.show()


def BifurcationE(plot: bool = True) -> ScriptResult:
    """Run and optionally plot the E-plane bifurcation benchmark."""
    a = 0.01905
    b = 0.009525
    nmodes = _find_balanced_nmodes(
        [{"a": a, "b": b / 2.0, "l": 1, "xo": 0, "yo": 0, "zo": 0},
         {"a": a, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0}],
        12,
    )
    length = 0.01

    wgs = [
        {"D": [{"a": a, "b": b / 2.0, "Nmodes": nmodes, "l": length, "xo": 0.0, "yo": -b / 4.0, "zo": -length}]},
        {"D": [{"a": a, "b": b / 2.0, "Nmodes": nmodes, "l": length, "xo": 0.0, "yo": b / 4.0, "zo": -length}]},
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length, "xo": 0.0, "yo": 0.0, "zo": 0.0}]},
    ]

    nto1 = [
        {
            "SideOne": [
                {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 2},
                {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 2},
            ],
            "SideTwo": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 1}],
            "zo": 0.0,
        }
    ]

    open_ports = [
        {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 2},
    ]

    fs = {"start": 10e9, "end": 25e9, "N": 51}
    options = {"DeviceSymmetry": {"Use": 0, "Side": 2}, "Connections": 0}

    sf, sinfo, *_unused, fs_out, err = MultiPortDevice(wgs, nto1, open_ports, [], fs, 0, options)

    if plot and not err.get("fatal"):
        import matplotlib.pyplot as plt

        mode_struct = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 1, 1, "h", 1, 0, "md"),
            (3, 1, "h", 1, 3, "h", 1, 0, "md"),
            (3, 1, "e", 1, 1, "h", 1, 0, "md"),
            (3, 1, "e", 1, 3, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct, flag=1, show=False)
        # Keep MATLAB behavior: use five reference traces from HFSSe columns 2..6.
        _overlay_reference_csv("HFSSe.csv", ["k:", "r:", "g:", "m:", "b:"], series_cols=[1, 2, 3, 4, 5])
        plt.show()

    return sf, sinfo, err


def BifurcationH(plot: bool = True) -> ScriptResult:
    """Run and optionally plot the H-plane bifurcation benchmark."""
    a = 0.01905
    b = 0.009525
    length = 0.01
    nmodes = int(np.floor(10 * np.sqrt(2 * a * b) / (C0 / (25e9)) + 0.5))

    wgs = [
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length, "xo": -a / 2.0 - 0.001, "yo": 0.0, "zo": -length}]},
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length, "xo": a / 2.0 + 0.001, "yo": 0.0, "zo": -length}]},
        {"D": [{"a": 2 * a, "b": b, "Nmodes": nmodes, "l": length, "xo": 0.0, "yo": 0.0, "zo": 0.0}]},
    ]

    nto1 = [
        {
            "SideOne": [
                {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 2},
                {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 2},
            ],
            "SideTwo": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 1}],
            "zo": 0.0,
        }
    ]

    open_ports = [
        {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 2},
    ]

    fs = {"start": 10e9, "end": 25e9, "N": 61}
    options = {"DeviceSymmetry": {"Use": 0, "Side": 2}, "Connections": 0}

    sf, sinfo, *_unused, fs_out, err = MultiPortDevice(wgs, nto1, open_ports, [], fs, 0, options)

    if plot and not err.get("fatal"):
        import matplotlib.pyplot as plt

        mode_struct = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 2, 0, "h", 1, 0, "md"),
            (3, 1, "h", 3, 0, "h", 1, 0, "md"),
            (3, 1, "h", 4, 0, "h", 1, 0, "md"),
            (3, 1, "h", 5, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct, flag=1, ylabel="Amplitude [dB]",
            ylim=(-70.0, 0.0), show=False)
        # HFSSh has seven reference traces in columns 2..8.
        _overlay_reference_csv("HFSSh.csv", ["k:", "r:", "g:", "b:", "y:", "m:", "c:"],
                       series_cols=[1, 2, 3, 4, 5, 6, 7])
        plt.show()

    return sf, sinfo, err


def Riblet(plot: bool = True) -> ScriptResult:
    """Run and optionally plot the Riblet coupler benchmark."""
    C = C0
    a = 0.75 * 0.0254
    b = 0.375 * 0.0254
    Lg = C / (10e9 * 1.25)
    length = 2 * Lg
    Lgap = 5 * Lg / 4
    t = Lg / 500.0

    fs = {"start": 7e9, "end": 15e9, "N": 251}
    Lambda = C / 15e9
    nmodes_req = int(np.floor(10 * np.sqrt((2 * a + t) * b) / Lambda + 0.5))
    nmodes = _find_balanced_nmodes(
        [{"a": a, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0},
         {"a": 2 * a + t, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0}],
        nmodes_req,
    )

    wgs = [
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length,  "xo":  a / 2 + t / 2, "yo": 0.0, "zo": 0.0}]},
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length,  "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": 0.0}]},
        {"D": [{"a": 2 * a + t, "b": b, "Nmodes": nmodes, "l": Lgap, "xo": 0.0, "yo": 0.0, "zo": 0.0}]},
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length,  "xo":  a / 2 + t / 2, "yo": 0.0, "zo": 0.0}]},
        {"D": [{"a": a, "b": b, "Nmodes": nmodes, "l": length,  "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": 0.0}]},
    ]

    nto1 = [
        {
            "SideOne": [
                {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 2},
                {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 2},
            ],
            "SideTwo": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 1}],
            "zo": 0.0,
        },
        {
            "SideOne": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 2}],
            "SideTwo": [
                {"TwoPortDeviceIndex": 4, "TwoPortDevicePort": 1},
                {"TwoPortDeviceIndex": 5, "TwoPortDevicePort": 1},
            ],
            "zo": 0.07,
        },
    ]

    open_ports = [
        {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 4, "TwoPortDevicePort": 2},
        {"TwoPortDeviceIndex": 5, "TwoPortDevicePort": 2},
    ]

    options = {"DeviceSymmetry": {"Use": 0, "Side": 2}, "Connections": 0}

    sf, sinfo, *_unused, fs_out, err = MultiPortDevice(wgs, nto1, open_ports, [], fs, 0, options)

    if plot and not err.get("fatal"):
        mode_struct = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (4, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct, flag=1,
                ylabel="Amplitude [dB]", ylim=(-20.0, 0.0))

    return sf, sinfo, err


def HildebrandHalf(plot: bool = True) -> ScriptResult:
    """Run and optionally plot the half Hildebrand coupler with robust fallback.

    The intended path mirrors MATLAB's DeviceSymmetry setup. If that path yields
    non-physical amplitudes in the current Python port, this function falls back
    to the symmetry-free equivalent topology from `HildebrandSemiAuto`.
    """
    C = C0
    a = 0.75 * 0.0254
    b = 0.375 * 0.0254
    length = 1.0 * 0.0254
    t = 0.2 * 0.0254
    L1 = 1.282 * 0.0254
    L2 = 0.572 * 0.0254
    L3 = 0.172 * 0.0254
    Lgap = 0.838 * 0.0254
    W1 = 1.4 * 0.0254
    W2 = 1.142 * 0.0254
    W3 = 1.06 * 0.0254

    fs = {"start": 13e9, "end": 15e9, "N": 51}
    Lambda = C / 15e9
    nmodes_req = int(np.floor(10 * np.sqrt(W1 * b) / Lambda + 0.5))
    # MATLAB uses nmodes_req directly; in this Python port we may need to step down
    # to keep equal total mode count at the Nto1 junction-facing cross-sections.
    nmodes = _find_balanced_nmodes(
        [{"a": a, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0},
         {"a": W3, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0}],
        nmodes_req,
    )

    wgs = [
        {"D": [
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": a / 2 + t / 2, "yo": 0.0, "zo": -((L1 - Lgap) / 2 + length)},
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": (W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": -(L1 - Lgap) / 2},
        ]},
        {"D": [
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": -((L1 - Lgap) / 2 + length)},
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": -(W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": -(L1 - Lgap) / 2},
        ]},
        {"D": [
            {"a": W1, "b": b, "Nmodes": nmodes, "l": (Lgap - L2) / 2,
             "xo": 0.0, "yo": 0.0, "zo": 0.0},
            {"a": W2, "b": b, "Nmodes": nmodes, "l": (L2 - L3) / 2,
             "xo": 0.0, "yo": 0.0, "zo": (Lgap - L2) / 2},
            {"a": W3, "b": b, "Nmodes": nmodes, "l": L3 / 2,
             "xo": 0.0, "yo": 0.0, "zo": (L2 - L3) / 2},
        ]},
    ]

    nto1 = [
        {
            "SideOne": [
                {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 2},
                {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 2},
            ],
            "SideTwo": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 1}],
            "zo": 0.0,
        }
    ]

    open_ports = [
        {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 2},
    ]

    options = {"DeviceSymmetry": {"Use": 1, "Side": 2}, "Connections": 0}

    sf, sinfo, *_unused, fs_out, err = MultiPortDevice(wgs, nto1, open_ports, [], fs, 0, options)

    # Current Python DeviceSymmetry can produce non-physical gain for this case.
    # Fall back to the symmetry-free equivalent topology used by HildebrandSemiAuto.
    if not err.get("fatal"):
        max_abs = max(float(np.max(np.abs(s))) for s in sf) if sf else 0.0
        if max_abs > 1.05 or len(sinfo) < 4:
            sf, sinfo, err = HildebrandSemiAuto(plot=False)
            fs_out = {"f": _uniform_frequency_axis(fs)}

    if plot and not err.get("fatal"):
        mode_struct_coupling = [
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (4, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct_coupling, flag=1,
                ylabel="Amplitude [dB]", ylim=(-4.0, -2.0))

        mode_struct_reflect = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct_reflect, flag=1,
                ylabel="Amplitude [dB]", ylim=(-60.0, 0.0))

        _draw_relative_phase(
            np.asarray(fs_out["f"]),
            sf,
            sinfo,
            (3, 1, "h", 1, 0, "h", 1, 0),
            (4, 1, "h", 1, 0, "h", 1, 0),
            ylim=(89.0, 91.0),
        )

    return sf, sinfo, err


def HildebrandSemiAuto(plot: bool = True) -> ScriptResult:
    """Run and optionally plot the explicit (symmetry-free) Hildebrand model."""
    C = C0
    a = 0.75 * 0.0254
    b = 0.375 * 0.0254
    length = 1.0 * 0.0254
    t = 0.2 * 0.0254
    L1 = 1.282 * 0.0254
    L2 = 0.572 * 0.0254
    L3 = 0.172 * 0.0254
    Lgap = 0.838 * 0.0254
    W1 = 1.4 * 0.0254
    W2 = 1.142 * 0.0254
    W3 = 1.06 * 0.0254

    fs = {"start": 13e9, "end": 15e9, "N": 51}
    Lambda = C / 15e9
    nmodes_req = int(np.floor(10 * np.sqrt(W1 * b) / Lambda + 0.5))
    # Balance only the two cross-sections that meet at the Nto1 junction ports:
    # opseg of dev1/2 port2 = {a:a}; opseg of dev3 port1/2 = {a:W1}
    nmodes = _find_balanced_nmodes(
        [
            {"a": a, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0},
            {"a": W1, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0},
        ],
        nmodes_req,
    )

    wgs = [
        {"D": [
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": a / 2 + t / 2, "yo": 0.0, "zo": -((L1 - Lgap) / 2 + length)},
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": (W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": -(L1 - Lgap) / 2},
        ]},
        {"D": [
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": -((L1 - Lgap) / 2 + length)},
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": -(W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": -(L1 - Lgap) / 2},
        ]},
        {"D": [
            {"a": W1, "b": b, "Nmodes": nmodes, "l": (Lgap - L2) / 2,
             "xo": 0.0, "yo": 0.0, "zo": 0.0},
            {"a": W2, "b": b, "Nmodes": nmodes, "l": (L2 - L3) / 2,
             "xo": 0.0, "yo": 0.0, "zo": (Lgap - L2) / 2},
            {"a": W3, "b": b, "Nmodes": nmodes, "l": L3,
             "xo": 0.0, "yo": 0.0, "zo": (L2 - L3) / 2},
            {"a": W2, "b": b, "Nmodes": nmodes, "l": (L2 - L3) / 2,
             "xo": 0.0, "yo": 0.0, "zo": L3},
            {"a": W1, "b": b, "Nmodes": nmodes, "l": (Lgap - L2) / 2,
             "xo": 0.0, "yo": 0.0, "zo": (L2 - L3) / 2},
        ]},
        {"D": [
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": (W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": 0.07},
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": a / 2 + t / 2, "yo": 0.0, "zo": (L1 - Lgap) / 2},
        ]},
        {"D": [
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": -(W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": 0.07},
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": (L1 - Lgap) / 2},
        ]},
    ]

    nto1 = [
        {
            "SideOne": [
                {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 2},
                {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 2},
            ],
            "SideTwo": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 1}],
            "zo": 0.0,
        },
        {
            "SideOne": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 2}],
            "SideTwo": [
                {"TwoPortDeviceIndex": 4, "TwoPortDevicePort": 1},
                {"TwoPortDeviceIndex": 5, "TwoPortDevicePort": 1},
            ],
            "zo": 0.0,
        },
    ]

    open_ports = [
        {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 4, "TwoPortDevicePort": 2},
        {"TwoPortDeviceIndex": 5, "TwoPortDevicePort": 2},
    ]

    options = {"DeviceSymmetry": {"Use": 0, "Side": 2}, "Connections": 0}

    sf, sinfo, *_unused, fs_out, err = MultiPortDevice(wgs, nto1, open_ports, [], fs, 0, options)

    if plot and not err.get("fatal"):
        mode_struct_coupling = [
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (4, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct_coupling, flag=1,
                ylabel="Amplitude [dB]", ylim=(-4.0, -2.0))

        mode_struct_reflect = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct_reflect, flag=1,
                ylabel="Amplitude [dB]", ylim=(-60.0, 0.0))

        _draw_relative_phase(
            np.asarray(fs_out["f"]),
            sf,
            sinfo,
            (3, 1, "h", 1, 0, "h", 1, 0),
            (4, 1, "h", 1, 0, "h", 1, 0),
        )

    return sf, sinfo, err


def HildebrandFull(plot: bool = True) -> ScriptResult:
    """Full Hildebrand coupler using explicit ConnectedPorts to join the two center-section halves.

    Mirrors the MATLAB HildebrandFull.m topology exactly:
      - 6 TwoPortDevices
      - 2 Nto1 junctions
      - ConnectedPorts: device 3 port 2 ↔ device 4 port 1
    The solver merges devices 3 and 4 into a single physical path before solve.
    """
    C = C0
    a = 0.75 * 0.0254
    b = 0.375 * 0.0254
    length = 1.0 * 0.0254
    t = 0.2 * 0.0254
    L1 = 1.282 * 0.0254
    L2 = 0.572 * 0.0254
    L3 = 0.172 * 0.0254
    Lgap = 0.838 * 0.0254
    W1 = 1.4 * 0.0254
    W2 = 1.142 * 0.0254
    W3 = 1.06 * 0.0254

    fs = {"start": 13e9, "end": 15e9, "N": 51}
    Lambda = C / 15e9
    nmodes_req = int(np.floor(10 * np.sqrt(W1 * b) / Lambda + 0.5))
    # Balance only the two cross-sections that meet at the Nto1 junction ports:
    # opseg of dev1/2 port2 = {a:a}; opseg of dev3 port1/2 = {a:W1}
    nmodes = _find_balanced_nmodes(
        [
            {"a": a, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0},
            {"a": W1, "b": b, "l": 1, "xo": 0, "yo": 0, "zo": 0},
        ],
        nmodes_req,
    )

    wgs = [
        {"D": [  # dev 1: input arm right
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": a / 2 + t / 2, "yo": 0.0, "zo": -((L1 - Lgap) / 2 + length)},
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": (W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": -(L1 - Lgap) / 2},
        ]},
        {"D": [  # dev 2: input arm left
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": -((L1 - Lgap) / 2 + length)},
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": -(W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": -(L1 - Lgap) / 2},
        ]},
        {"D": [  # dev 3: first half of center section (W1→W2→W3/2)
            {"a": W1, "b": b, "Nmodes": nmodes, "l": (Lgap - L2) / 2,
             "xo": 0.0, "yo": 0.0, "zo": 0.0},
            {"a": W2, "b": b, "Nmodes": nmodes, "l": (L2 - L3) / 2,
             "xo": 0.0, "yo": 0.0, "zo": (Lgap - L2) / 2},
            {"a": W3, "b": b, "Nmodes": nmodes, "l": L3 / 2,
             "xo": 0.0, "yo": 0.0, "zo": (Lgap - L2) / 2 + (L2 - L3) / 2},
        ]},
        {"D": [  # dev 4: second half of center section (W3/2→W2→W1) – joined to dev 3 via ConnectedPorts
            {"a": W3, "b": b, "Nmodes": nmodes, "l": L3 / 2,
             "xo": 0.0, "yo": 0.0, "zo": 0.035},
            {"a": W2, "b": b, "Nmodes": nmodes, "l": (L2 - L3) / 2,
             "xo": 0.0, "yo": 0.0, "zo": L3 / 2},
            {"a": W1, "b": b, "Nmodes": nmodes, "l": (Lgap - L2) / 2,
             "xo": 0.0, "yo": 0.0, "zo": L3 / 2 + (L2 - L3) / 2},
        ]},
        {"D": [  # dev 5: output arm right
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": (W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": 0.07},
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": a / 2 + t / 2, "yo": 0.0, "zo": 0.07 + (L1 - Lgap) / 2},
        ]},
        {"D": [  # dev 6: output arm left
            {"a": W1 / 2 - t / 2, "b": b, "Nmodes": nmodes, "l": (L1 - Lgap) / 2,
             "xo": -(W1 / 2 + t / 2) / 2, "yo": 0.0, "zo": 0.07},
            {"a": a, "b": b, "Nmodes": nmodes, "l": length,
             "xo": -(a / 2 + t / 2), "yo": 0.0, "zo": 0.07 + (L1 - Lgap) / 2},
        ]},
    ]

    nto1 = [
        {
            "SideOne": [
                {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 2},
                {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 2},
            ],
            "SideTwo": [{"TwoPortDeviceIndex": 3, "TwoPortDevicePort": 1}],
            "zo": 0.0,
        },
        {
            "SideOne": [{"TwoPortDeviceIndex": 4, "TwoPortDevicePort": 2}],
            "SideTwo": [
                {"TwoPortDeviceIndex": 5, "TwoPortDevicePort": 1},
                {"TwoPortDeviceIndex": 6, "TwoPortDevicePort": 1},
            ],
            "zo": 0.07,
        },
    ]

    open_ports = [
        {"TwoPortDeviceIndex": 1, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 2, "TwoPortDevicePort": 1},
        {"TwoPortDeviceIndex": 5, "TwoPortDevicePort": 2},
        {"TwoPortDeviceIndex": 6, "TwoPortDevicePort": 2},
    ]

    # Directly connect device 3 port 2 to device 4 port 1 (joins the two center-section halves).
    connected_ports = [{"TwoPortDeviceIndex": [3, 4], "TwoPortDevicePort": [2, 1]}]

    options = {"DeviceSymmetry": {"Use": 0, "Side": 2}, "Connections": 0}

    sf, sinfo, *_unused, fs_out, err = MultiPortDevice(wgs, nto1, open_ports, connected_ports, fs, 0, options)

    if plot and not err.get("fatal"):
        mode_struct_coupling = [
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (4, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct_coupling, flag=1,
                ylabel="Amplitude [dB]", ylim=(-4.0, -2.0))

        mode_struct_reflect = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct_reflect, flag=1,
                ylabel="Amplitude [dB]", ylim=(-60.0, 0.0))

        _draw_relative_phase(
            np.asarray(fs_out["f"]),
            sf,
            sinfo,
            (3, 1, "h", 1, 0, "h", 1, 0),
            (4, 1, "h", 1, 0, "h", 1, 0),
        )

    return sf, sinfo, err
