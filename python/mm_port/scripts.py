from __future__ import annotations

from typing import Any

import numpy as np

from .core import C0, GSMDraw, MultiPortDevice, OrderModes


def _find_balanced_nmodes(segments_specs: list[dict], nmodes_max: int) -> int:
    """Return the largest n ≤ nmodes_max where all segments yield the same total mode count."""
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


def BifurcationE(plot: bool = True) -> tuple[list[np.ndarray], list[dict[str, Any]], dict[str, Any]]:
    a = 0.01905
    b = 0.009525
    nmodes = 10
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
        mode_struct = [
            (1, 1, "h", 1, 0, "h", 1, 0, "md"),
            (2, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 1, 0, "h", 1, 0, "md"),
            (3, 1, "h", 1, 1, "h", 1, 0, "md"),
            (3, 1, "e", 1, 1, "h", 1, 0, "md"),
        ]
        GSMDraw(np.asarray(fs_out["f"]), sf, sinfo, mode_struct, flag=1)

    return sf, sinfo, err


def BifurcationH(plot: bool = True) -> tuple[list[np.ndarray], list[dict[str, Any]], dict[str, Any]]:
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
                ylim=(-70.0, 0.0))

    return sf, sinfo, err


def Riblet(plot: bool = True) -> tuple[list[np.ndarray], list[dict[str, Any]], dict[str, Any]]:
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


def HildebrandHalf(plot: bool = True) -> tuple[list[np.ndarray], list[dict[str, Any]], dict[str, Any]]:
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

    return sf, sinfo, err


def HildebrandSemiAuto(*args: Any, **kwargs: Any) -> Any:
    raise NotImplementedError("HildebrandSemiAuto pure-Python script port is not implemented yet")


def HildebrandFull(*args: Any, **kwargs: Any) -> Any:
    raise NotImplementedError("HildebrandFull pure-Python script port is not implemented yet")
