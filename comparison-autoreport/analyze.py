#!/usr/bin/env python3
"""Build outliers, plots, and summary.json from a FORTRAN/Cpp dump pair.

The dump HEADER/FOOTER (kinematics, m/m2, call counts) is the source of
truth for the report. CLI flags fill in anything the dump does not carry
(seed, stream sha256, toolchain).
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "harness" / "python"))

from compare import parse_events, parse_footer, parse_header, rel_delta  # noqa: E402
from plot_deltas import notable, plot_deltas  # noqa: E402

from lib import (  # noqa: E402
    abs3,
    dump_diff,
    dump_json,
    ecm_pcm,
    fmt,
    hex32,
    k2_from,
    load_json,
    load_quads,
    mag3,
    parse_meta,
    write_outliers,
)

MEV = 1.0e3
M_DEFAULT = 0.511000e-3
M2_DEFAULT = 0.261112e-6
ELAB_DEFAULT = 45.0
PHYS_REL = 1.0e-6


def linspace(a: float, b: float, nbin: int) -> list[float]:
    if b == a:
        b = a + 1e-12
    w = (b - a) / nbin
    return [a + i * w for i in range(nbin + 1)]


def hist(xs: list[float], edges: list[float]) -> list[int]:
    counts = [0] * (len(edges) - 1)
    for x in xs:
        if x <= edges[0]:
            counts[0] += 1
            continue
        if x >= edges[-1]:
            counts[-1] += 1
            continue
        lo, hi = 0, len(edges) - 1
        while lo < hi - 1:
            mid = (lo + hi) // 2
            if x < edges[mid]:
                hi = mid
            else:
                lo = mid
        counts[lo] += 1
    return counts


def overlay_figure(
    out_png: Path,
    eg_f: list[float],
    eg_c: list[float],
    k2_f: list[float],
    k2_c: list[float],
    n: int,
    n_call_last: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eg_all = eg_f + eg_c
    k2_all = k2_f + k2_c
    eg_edges = linspace(min(eg_all), max(eg_all) * 1.001, 50)
    k2_edges = linspace(min(k2_all), max(k2_all) * 1.001, 50)
    w_eg = (eg_edges[1] - eg_edges[0]) * MEV
    w_k2 = (k2_edges[1] - k2_edges[0]) * MEV
    cx_eg = [0.5 * (eg_edges[i] + eg_edges[i + 1]) * MEV
             for i in range(len(eg_edges) - 1)]
    cx_k2 = [0.5 * (k2_edges[i] + k2_edges[i + 1]) * MEV
             for i in range(len(k2_edges) - 1)]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), layout="constrained")
    ax = axes[0]
    ax.bar(cx_eg, hist(eg_f, eg_edges), width=w_eg, align="center", fill=False,
           edgecolor="C0", linewidth=1.4, label="FORTRAN")
    ax.bar(cx_eg, hist(eg_c, eg_edges), width=w_eg, align="center", fill=False,
           edgecolor="C3", linewidth=1.0, linestyle="--", label="Cpp")
    ax.set_title("Radiated photon energy (CM)")
    ax.set_xlabel(r"$E_\gamma$ (MeV)")
    ax.set_ylabel("Radiative events")
    ax.legend(frameon=False)

    ax = axes[1]
    ax.bar(cx_k2, hist(k2_f, k2_edges), width=w_k2, align="center", fill=False,
           edgecolor="C0", linewidth=1.4, label="FORTRAN")
    ax.bar(cx_k2, hist(k2_c, k2_edges), width=w_k2, align="center", fill=False,
           edgecolor="C3", linewidth=1.0, linestyle="--", label="Cpp")
    ax.set_title("Outgoing electron $|k_2|$ (CM)")
    ax.set_xlabel(r"$|k_2|$ (MeV)")
    ax.set_ylabel("Radiative events")
    ax.legend(frameon=False)

    fig.suptitle(
        f"MERADGEN FORTRAN vs Cpp  ·  {n:,} radiative events  ·  "
        f"shared stream through call {n_call_last}",
        fontsize=10,
    )
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_bin_tsv(path: Path, centers: list[float], a: list[int], b: list[int],
                  unit: str) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"center_{unit}", "fortran", "cpp"])
        for c, fa, ca in zip(centers, a, b):
            w.writerow([f"{c:.8e}", fa, ca])


def first_quad_preview(quads: list[list[float]]) -> list[float]:
    if not quads:
        return []
    return [round(x, 4) for x in quads[0]]


def pick_kin(header: dict, key: str, fallback: float) -> float:
    val = header.get(key)
    return float(val) if val is not None else fallback


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Analyze a FORTRAN vs Cpp dump pair into report inputs."
    )
    ap.add_argument("--run-dir", type=Path, required=True,
                    help="directory holding dumps, quads, and config.json")
    ap.add_argument("--fortran", type=Path)
    ap.add_argument("--cpp", type=Path)
    ap.add_argument("--quads", type=Path)
    ap.add_argument("--elab", type=float, default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    out_dir = args.run_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    config_path = out_dir / "config.json"
    config = load_json(config_path) if config_path.is_file() else {}

    out_f = args.fortran or Path(config.get("fortran_dump", out_dir / "fortran_radiative.txt"))
    out_c = args.cpp or Path(config.get("cpp_dump", out_dir / "cpp_radiative.txt"))
    quads_path = args.quads or Path(config.get("quads", out_dir / "quads.txt"))
    if not out_f.is_file() or not out_c.is_file():
        print(f"FAIL missing dumps:\n  {out_f}\n  {out_c}", file=sys.stderr)
        return 1

    hf = parse_header(out_f)
    hc = parse_header(out_c)
    ff = parse_footer(out_f)
    fc = parse_footer(out_c)
    header = dict(hf)
    if hf and hc:
        for key in ("elab", "thetacm", "phi", "pl"):
            if key in hf and key in hc and hf[key] != hc[key]:
                print(f"FAIL dump HEADER {key} mismatch "
                      f"fortran={hf[key]} cpp={hc[key]}", file=sys.stderr)
                return 1

    elab = args.elab
    if elab is None:
        elab = pick_kin(header, "elab", float(config.get("elab", ELAB_DEFAULT)))
    m = pick_kin(header, "m", M_DEFAULT)
    m2 = pick_kin(header, "m2", M2_DEFAULT)
    if "ecm" in header and "pcm" in header:
        ecm, pcm = float(header["ecm"]), float(header["pcm"])
    else:
        ecm, pcm = ecm_pcm(elab, m, m2)

    fe = parse_events(out_f)
    ce = parse_events(out_c)
    if len(fe) != len(ce):
        print(f"FAIL radiative count fortran={len(fe)} cpp={len(ce)}", file=sys.stderr)
        return 1
    if not fe:
        print("FAIL no events in dumps", file=sys.stderr)
        return 1

    quads = load_quads(quads_path if quads_path.is_file() else None)
    meta = parse_meta(quads_path.with_suffix(".meta")) if quads_path.is_file() else {}

    events: list[int] = []
    eg_f: list[float] = []
    eg_c: list[float] = []
    k2m_f: list[float] = []
    k2m_c: list[float] = []
    dE_g: list[float] = []
    dp_g: list[float] = []
    dE_e: list[float] = []
    dp_e: list[float] = []
    dmag_e: list[float] = []
    outliers: list[dict] = []
    n_ich = 0
    n_rel_1e6 = 0
    n_bit_phirad = 0
    n_bit_vprad = 0
    worst_g = (0.0, 0)
    worst_e = (0.0, 0)

    for f, c in zip(fe, ce):
        if f["index"] != c["index"]:
            print(f"FAIL index mismatch {f['index']} vs {c['index']}", file=sys.stderr)
            return 1
        ich_f = f["kin"]["ich"] if f.get("kin") else 1
        ich_c = c["kin"]["ich"] if c.get("kin") else 1
        if ich_f != ich_c:
            n_ich += 1
        pg_f, pg_c = f["phirad"], c["phirad"]
        pe_f = k2_from(f["vprad"], f["phirad"], ecm, pcm)
        pe_c = k2_from(c["vprad"], c["phirad"], ecm, pcm)
        de_g = pg_c[3] - pg_f[3]
        d3_g = abs3(pg_c, pg_f)
        de_e = pe_c[3] - pe_f[3]
        d3_e = abs3(pe_c, pe_f)
        dabs_e = mag3(pe_c) - mag3(pe_f)
        rg = rel_delta(pg_f[3], pg_c[3])
        re_ = rel_delta(pe_f[3], pe_c[3])
        if any(rel_delta(a, b) > PHYS_REL for a, b in
               zip(f["vprad"] + f["phirad"], c["vprad"] + c["phirad"])):
            n_rel_1e6 += 1
        if f["phirad"] != c["phirad"]:
            n_bit_phirad += 1
        if f["vprad"] != c["vprad"]:
            n_bit_vprad += 1

        events.append(f["index"])
        eg_f.append(pg_f[3])
        eg_c.append(pg_c[3])
        k2m_f.append(mag3(pe_f))
        k2m_c.append(mag3(pe_c))
        dE_g.append(de_g)
        dp_g.append(d3_g)
        dE_e.append(de_e)
        dp_e.append(d3_e)
        dmag_e.append(dabs_e)

        mag_g = max(abs(de_g), d3_g)
        mag_e = max(abs(de_e), d3_e, abs(dabs_e))
        if mag_g > worst_g[0]:
            worst_g = (mag_g, f["index"])
        if mag_e > worst_e[0]:
            worst_e = (mag_e, f["index"])

        keep = ich_f != ich_c
        keep = keep or notable(abs(de_g), rg, 0.0, 0.0)
        keep = keep or notable(d3_g, d3_g / max(mag3(pg_f), mag3(pg_c), 1e-30),
                               0.0, 0.0)
        keep = keep or notable(abs(de_e), re_, 0.0, 0.0)
        keep = keep or notable(d3_e, d3_e / max(mag3(pe_f), mag3(pe_c), 1e-30),
                               0.0, 0.0)
        keep = keep or notable(abs(dabs_e),
                               abs(dabs_e) / max(mag3(pe_f), mag3(pe_c), 1e-30),
                               0.0, 0.0)
        if not keep:
            continue
        if quads and 1 <= f["index"] <= len(quads):
            rnd = quads[f["index"] - 1]
        elif f.get("random"):
            rnd = f["random"]
        else:
            rnd = [float("nan")] * 4
        outliers.append({
            "event": f["index"],
            "ich_f": ich_f,
            "ich_c": ich_c,
            "r0": rnd[0], "r1": rnd[1], "r2": rnd[2], "r3": rnd[3],
            "vp_fx": f["vprad"][0], "vp_fy": f["vprad"][1],
            "vp_fz": f["vprad"][2], "vp_fE": f["vprad"][3],
            "vp_cx": c["vprad"][0], "vp_cy": c["vprad"][1],
            "vp_cz": c["vprad"][2], "vp_cE": c["vprad"][3],
            "pg_fx": pg_f[0], "pg_fy": pg_f[1], "pg_fz": pg_f[2], "pg_fE": pg_f[3],
            "pg_cx": pg_c[0], "pg_cy": pg_c[1], "pg_cz": pg_c[2], "pg_cE": pg_c[3],
            "pe_fx": pe_f[0], "pe_fy": pe_f[1], "pe_fz": pe_f[2], "pe_fE": pe_f[3],
            "pe_cx": pe_c[0], "pe_cy": pe_c[1], "pe_cz": pe_c[2], "pe_cE": pe_c[3],
            "dE_g": de_g, "dp_g": d3_g, "dE_e": de_e, "dp_e": d3_e,
            "dabs_k2": dabs_e, "rel_Eg": rg, "rel_Ee": re_,
        })

    outliers.sort(key=lambda r: max(abs(r["dE_g"]), r["dp_g"], abs(r["dE_e"]),
                                    r["dp_e"], abs(r["dabs_k2"])), reverse=True)
    n = len(fe)
    n_zero_g = sum(1 for a, b in zip(dE_g, dp_g) if a == 0.0 and b == 0.0)
    n_zero_e = sum(1 for a, b, c in zip(dE_e, dp_e, dmag_e)
                   if a == 0.0 and b == 0.0 and c == 0.0)
    nonzero_g = [i for i, (a, b) in enumerate(zip(dE_g, dp_g)) if a != 0.0 or b != 0.0]
    nonzero_e = [i for i, (a, b, c) in enumerate(zip(dE_e, dp_e, dmag_e))
                 if a != 0.0 or b != 0.0 or c != 0.0]
    n_call_last = events[-1] if events else 0
    n_calls = int(ff.get("calls", fc.get("calls", n_call_last)))
    n_rad_footer = int(ff.get("radiative", n))
    if ff and fc and ff.get("calls") != fc.get("calls"):
        print(f"FAIL call-count mismatch fortran={ff.get('calls')} "
              f"cpp={fc.get('calls')}", file=sys.stderr)
        return 1

    tsv = out_dir / "outliers.tsv"
    write_outliers(tsv, outliers)

    eg_edges = linspace(min(eg_f + eg_c), max(eg_f + eg_c) * 1.001, 50)
    k2_edges = linspace(min(k2m_f + k2m_c), max(k2m_f + k2m_c) * 1.001, 50)
    eg_hf, eg_hc = hist(eg_f, eg_edges), hist(eg_c, eg_edges)
    k2_hf, k2_hc = hist(k2m_f, k2_edges), hist(k2m_c, k2_edges)
    cx_eg = [0.5 * (eg_edges[i] + eg_edges[i + 1]) * MEV for i in range(50)]
    cx_k2 = [0.5 * (k2_edges[i] + k2_edges[i + 1]) * MEV for i in range(50)]
    write_bin_tsv(out_dir / "bins_egamma.tsv", cx_eg, eg_hf, eg_hc, "MeV")
    write_bin_tsv(out_dir / "bins_k2.tsv", cx_k2, k2_hf, k2_hc, "MeV")

    n1 = sum(1 for r in outliers if r["dE_g"] != 0.0)
    n2 = sum(1 for r in outliers if r["dE_g"] != 0.0 or r["dp_g"] != 0.0)
    n3 = sum(1 for r in outliers if r["dE_e"] != 0.0)
    n4 = sum(1 for r in outliers if r["dE_e"] != 0.0 or r["dp_e"] != 0.0
             or r["dabs_k2"] != 0.0)

    hex_example = None
    if outliers:
        r0 = outliers[0]
        label, a, b, ulps, ad = dump_diff(r0)
        rel = rel_delta(a, b) if a != 0.0 or b != 0.0 else 0.0
        hex_example = {
            "event": r0["event"],
            "component": label,
            "fortran_hex": hex32(a),
            "cpp_hex": hex32(b),
            "abs_delta_GeV": ad,
            "ulps": ulps,
            "rel": rel,
        }

    mismatches = []
    for r in outliers:
        label, a, b, ulps, ad = dump_diff(r)
        mismatches.append({
            "event": r["event"],
            "component": label,
            "fortran_hex": hex32(a),
            "cpp_hex": hex32(b),
            "abs_delta_GeV": ad,
            "ulps": ulps,
        })

    toolchain = ""
    tpath = out_dir / "toolchain.txt"
    if tpath.is_file():
        toolchain = tpath.read_text()

    n_quads = int(meta.get("n", config.get("n_quads", len(quads))))
    seed = meta.get("seed", config.get("seed"))
    frac = (n_rad_footer / n_calls) if n_calls else 0.0

    summary = {
        "date": str(config.get("date", date.today().isoformat())),
        "command": config.get("command", []),
        "seed": seed,
        "n_quads": n_quads,
        "n_radiative_requested": config.get("n_radiative"),
        "n_radiative": n,
        "n_calls": n_calls,
        "last_stream_index": n_call_last,
        "radiative_fraction": frac,
        "quads_sha256": meta.get("sha256"),
        "first_quad": quads[0] if quads else None,
        "first_quad_preview": first_quad_preview(quads),
        "kinematics": {
            "elab": elab,
            "thetacm": pick_kin(header, "thetacm", float(config.get("thetacm", 90.0))),
            "phi": pick_kin(header, "phi", float(config.get("phi", 10.0))),
            "pl": pick_kin(header, "pl", float(config.get("pl", -1.0))),
            "m": m,
            "m2": m2,
            "ecm": ecm,
            "pcm": pcm,
            "vpgen": header.get("vpgen"),
        },
        "dump_header_fortran": hf,
        "dump_header_cpp": hc,
        "dump_footer_fortran": ff,
        "dump_footer_cpp": fc,
        "prec": header.get("mode", config.get("prec", "full")),
        "dump_policy": header.get("dump_policy", "radiative_only"),
        "toolchain_text": toolchain,
        "flags": "-O2 -ffp-contract=off -fno-fast-math",
        "ich_mismatch": n_ich,
        "vprad_any_diff": n_bit_vprad,
        "phirad_any_diff": n_bit_phirad,
        "events_above_rel_1e-6": n_rel_1e6,
        "physics_rel_bar": PHYS_REL,
        "photon_exact_zero": n_zero_g,
        "electron_exact_zero": n_zero_e,
        "photon_nonzero": len(nonzero_g),
        "electron_nonzero": len(nonzero_e),
        "n_dump_mismatch": len(outliers),
        "panel1_dEg": n1,
        "panel2_photon": n2,
        "panel3_dEe": n3,
        "panel4_electron": n4,
        "overlay_bins_identical_egamma": eg_hf == eg_hc,
        "overlay_bins_identical_k2": k2_hf == k2_hc,
        "worst_photon": {"abs_GeV": worst_g[0], "event": worst_g[1]},
        "worst_electron": {"abs_GeV": worst_e[0], "event": worst_e[1]},
        "hex_example": hex_example,
        "mismatches": mismatches,
        "files": {
            "fortran_dump": str(out_f),
            "cpp_dump": str(out_c),
            "quads": str(quads_path) if quads_path else None,
            "outliers_tsv": str(tsv),
            "overlay_png": str(out_dir / "overlay_radiative.png"),
            "delta_png": str(out_dir / "delta_photon_electron.png"),
        },
        "trees": {
            "fortran": "meradgen-fortran/",
            "cpp": "meradgen-cpp/",
            "cpp_final": "meradgen-cpp-final/",
        },
    }
    dump_json(out_dir / "summary.json", summary)

    lines = [
        f"radiative_events={n}",
        f"last_stream_index={n_call_last}",
        f"n_calls={n_calls}",
        f"radiative_fraction={frac:.6f}",
        f"ich_mismatch={n_ich}",
        f"vprad_any_diff={n_bit_vprad}",
        f"phirad_any_diff={n_bit_phirad}",
        f"events_above_rel_1e-6={n_rel_1e6}",
        f"photon_exact_zero={n_zero_g}",
        f"electron_exact_zero={n_zero_e}",
        f"photon_nonzero={len(nonzero_g)}",
        f"electron_nonzero={len(nonzero_e)}",
        f"outliers={len(outliers)}",
        f"elab={elab}",
        f"thetacm={summary['kinematics']['thetacm']}",
        f"phi={summary['kinematics']['phi']}",
        f"pl={summary['kinematics']['pl']}",
        f"overlay_bins_identical_egamma={eg_hf == eg_hc}",
        f"overlay_bins_identical_k2={k2_hf == k2_hc}",
        f"summary_json={out_dir / 'summary.json'}",
    ]
    (out_dir / "summary.txt").write_text("\n".join(lines) + "\n")
    for line in lines:
        print(line)

    if not args.no_plot:
        overlay_figure(
            out_dir / "overlay_radiative.png",
            eg_f, eg_c, k2m_f, k2m_c, n, n_call_last,
        )
        print(f"plot={out_dir / 'overlay_radiative.png'}")
        plot_deltas(
            out_dir / "delta_photon_electron.png",
            dE_g, dp_g, dE_e, dp_e, dmag_e, events, nonzero_g, nonzero_e,
            n, n_zero_g, n_zero_e,
        )
        print(f"plot={out_dir / 'delta_photon_electron.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
