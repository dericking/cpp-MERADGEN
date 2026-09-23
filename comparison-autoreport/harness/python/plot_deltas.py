#!/usr/bin/env python3
"""Plot and dump FORTRAN vs C++ deltas for the real photon and outgoing electron.

Both generators consume one shared random stream (generate_stream.py / quads.txt).
Full event dumps stay in _scratch/; this writes only a summary plus the events
whose photon or electron four-vectors differ enough to investigate.

Critical path is stdlib. matplotlib is optional and used only for the figures.
"""
from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import parse_events, rel_delta  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
HARNESS = Path(__file__).resolve().parents[1]
VALROOT = Path(__file__).resolve().parents[2]

# Parity-tree literals (meradgen-fortran/run.f DATA; meradgen-cpp/globals.cpp).
# Used only to reconstruct k2 = k1+p1-p2-γ from VPRAD/PHIRAD. Must match the
# driver that produced the dumps.
M = 0.511000e-3
M2 = 0.261112e-6
ELAB_DEFAULT = 45.0


def ecm_pcm(elab: float) -> tuple[float, float]:
    ecm = math.sqrt(2.0 * M * (elab + M)) / 2.0
    pcm = math.sqrt(ecm * ecm - M2)
    return ecm, pcm


def k2_from(vprad: list[float], phirad: list[float], ecm: float, pcm: float) -> list[float]:
    """Outgoing beam electron (px, py, pz, E) in the CM frame."""
    return [
        -vprad[0] - phirad[0],
        -vprad[1] - phirad[1],
        pcm - vprad[2] - phirad[2],
        ecm - vprad[3] - phirad[3],
    ]


def abs3(a: list[float], b: list[float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def mag3(v: list[float]) -> float:
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def fmt(x: float) -> str:
    return f"{x:.16e}"


def load_quads(path: Path | None, n: int) -> list[list[float] | None]:
    rows: list[list[float] | None] = [None] * n
    if path is None or not path.is_file():
        return rows
    with path.open() as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            parts = line.split()
            if len(parts) >= 4:
                rows[i] = [float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])]
    return rows


def notable(absd: float, reld: float, keep_abs: float, keep_rel: float) -> bool:
    if absd == 0.0:
        return False
    if keep_abs == 0.0 and keep_rel == 0.0:
        return True
    return absd >= keep_abs or reld >= keep_rel


def run_generators(campaign: Path, seed: int, n: int, quads: Path | None) -> Path:
    cmd = [
        sys.executable,
        str(HARNESS / "python" / "run.py"),
        "--campaign",
        campaign.name,
        "--prec",
        "full",
    ]
    if seed is not None and n is not None:
        cmd.extend(["--seed", str(seed), "--n", str(n)])
    elif quads is not None:
        cmd.extend(["--quads", str(quads)])
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=False)
    return campaign / "_scratch"


def write_summary(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def hist_range(xs: list[float]) -> tuple[float, float]:
    lo, hi = min(xs), max(xs)
    if lo == hi:
        pad = 1e-15 if lo == 0.0 else abs(lo) * 1e-3 + 1e-18
        return lo - pad, hi + pad
    pad = 0.08 * (hi - lo)
    return lo - pad, hi + pad


def plot_deltas(
    out_png: Path,
    dE_g: list[float],
    dp_g: list[float],
    dE_e: list[float],
    dp_e: list[float],
    dmag_e: list[float],
    events: list[int],
    nonzero_g: list[int],
    nonzero_e: list[int],
    n: int,
    n_zero_g: int,
    n_zero_e: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.6), layout="constrained")

    xmax = max(events) if events else n

    ax = axes[0, 0]
    lo, hi = hist_range(dE_g)
    ax.hist(dE_g, bins=80, range=(lo, hi), histtype="step", color="C0", linewidth=1.2)
    ax.set_yscale("log")
    ax.set_title(r"(1) Photon $\Delta E_\gamma$ (Cpp $-$ FORTRAN)")
    ax.set_xlabel(r"$\Delta E_\gamma$ (GeV)")
    ax.set_ylabel("Events")
    n_de_g0 = sum(1 for x in dE_g if x == 0.0)
    ax.text(
        0.02,
        0.97,
        f"{n_de_g0}/{n} with $\\Delta E=0$\n{n - n_de_g0} nonzero (Table 1)",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color="0.35",
    )

    ax = axes[0, 1]
    y_g = [max(abs(dE_g[i]), dp_g[i]) for i in nonzero_g]
    if nonzero_g:
        ax.plot(
            [events[i] for i in nonzero_g],
            y_g,
            ".",
            ms=5,
            color="C0",
        )
        ax.set_yscale("log")
    else:
        ax.text(0.5, 0.5, "no nonzero photon deltas", ha="center", va="center",
                transform=ax.transAxes, color="0.4")
    ax.set_xlim(1, xmax)
    ax.set_title(r"(2) Photon max($|\Delta E_\gamma|,|\Delta\vec{p}_\gamma|$) vs event")
    ax.set_xlabel("Event index (1-based)")
    ax.set_ylabel(r"max $|\Delta|$ (GeV)")

    ax = axes[1, 0]
    lo, hi = hist_range(dE_e)
    ax.hist(dE_e, bins=80, range=(lo, hi), histtype="step", color="C3", linewidth=1.2)
    ax.set_yscale("log")
    ax.set_title(r"(3) Electron $\Delta E_{k_2}$ (Cpp $-$ FORTRAN)")
    ax.set_xlabel(r"$\Delta E_{k_2}$ (GeV)")
    ax.set_ylabel("Events")
    n_de_e0 = sum(1 for x in dE_e if x == 0.0)
    ax.text(
        0.02,
        0.97,
        f"{n_de_e0}/{n} with $\\Delta E=0$\n{n - n_de_e0} nonzero (Table 3)",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color="0.35",
    )

    ax = axes[1, 1]
    y_e = [max(abs(dE_e[i]), dp_e[i], abs(dmag_e[i])) for i in nonzero_e]
    if nonzero_e:
        ax.plot(
            [events[i] for i in nonzero_e],
            y_e,
            ".",
            ms=5,
            color="C3",
        )
        ax.set_yscale("log")
    else:
        ax.text(0.5, 0.5, "no nonzero electron deltas", ha="center", va="center",
                transform=ax.transAxes, color="0.4")
    ax.set_xlim(1, xmax)
    ax.set_title(r"(4) Electron max($|\Delta E|,|\Delta\vec{k}_2|,|\Delta|k_2||$) vs event")
    ax.set_xlabel("Event index (1-based)")
    ax.set_ylabel(r"max $|\Delta|$ (GeV)")

    fig.suptitle("MERADGEN Cpp $-$ FORTRAN  ·  same random stream", fontsize=11)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Photon and outgoing-electron deltas, FORTRAN vs C++, shared seed."
    )
    ap.add_argument("--campaign", help="dated dir under validation_checks_new/")
    ap.add_argument("--fortran", type=Path, help="fortran_output.txt (overrides campaign scratch)")
    ap.add_argument("--cpp", type=Path, help="cpp_output.txt (overrides campaign scratch)")
    ap.add_argument("--quads", type=Path, help="shared quads.txt (full-precision rand4 for outliers)")
    ap.add_argument("--out-dir", type=Path, help="where to write TSV/PNG (default: campaign/_scratch)")
    ap.add_argument("--seed", type=int, help="with --run --n, generate the shared stream")
    ap.add_argument("--n", type=int, help="event count when generating a stream")
    ap.add_argument("--run", action="store_true",
                    help="build and run both generators on the shared stream first")
    ap.add_argument("--elab", type=float, default=ELAB_DEFAULT,
                    help="lab energy used to reconstruct k2 (must match the driver)")
    ap.add_argument("--keep-abs", type=float, default=0.0,
                    help="keep outlier if any |ΔE| or |Δp| is ≥ this (GeV). "
                         "0 keeps every nonzero delta")
    ap.add_argument("--keep-rel", type=float, default=0.0,
                    help="also keep if relative |ΔE| is ≥ this")
    ap.add_argument("--no-plot", action="store_true", help="write data files only")
    args = ap.parse_args()

    campaign = VALROOT / args.campaign if args.campaign else None
    if args.run:
        if campaign is None:
            raise SystemExit("--run requires --campaign")
        if args.seed is None and args.quads is None and args.n is None:
            raise SystemExit("--run needs --seed/--n or --quads")
        campaign.mkdir(parents=True, exist_ok=True)
        scratch = run_generators(campaign, args.seed, args.n, args.quads)
    else:
        scratch = campaign / "_scratch" if campaign is not None else None

    out_f = args.fortran
    out_c = args.cpp
    if out_f is None or out_c is None:
        if scratch is None:
            raise SystemExit("need --fortran/--cpp or --campaign")
        out_f = out_f or (scratch / "fortran_output.txt")
        out_c = out_c or (scratch / "cpp_output.txt")
    if not out_f.is_file() or not out_c.is_file():
        raise SystemExit(f"missing dumps:\n  {out_f}\n  {out_c}")

    out_dir = args.out_dir
    if out_dir is None:
        out_dir = scratch if scratch is not None else out_f.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    quads_path = args.quads
    if quads_path is None and scratch is not None:
        cand = scratch / "quads.txt"
        if cand.is_file():
            quads_path = cand

    fe = parse_events(out_f)
    ce = parse_events(out_c)
    if len(fe) != len(ce):
        print(f"FAIL event count fortran={len(fe)} cpp={len(ce)}", file=sys.stderr)
        return 1

    ecm, pcm = ecm_pcm(args.elab)
    quads = load_quads(quads_path, len(fe))

    events: list[int] = []
    dE_g: list[float] = []
    dp_g: list[float] = []
    dE_e: list[float] = []
    dp_e: list[float] = []
    dmag_e: list[float] = []
    rel_g: list[float] = []
    rel_e: list[float] = []
    outliers: list[dict] = []
    n_rad = 0
    n_ich = 0
    worst_g = (0.0, 0, "none")
    worst_e = (0.0, 0, "none")

    for f, c in zip(fe, ce):
        if f["index"] != c["index"]:
            print(f"FAIL index mismatch {f['index']} vs {c['index']}", file=sys.stderr)
            return 1
        if f["vprad"] is None or c["vprad"] is None or f["phirad"] is None or c["phirad"] is None:
            print(f"FAIL missing vectors event={f['index']}", file=sys.stderr)
            return 1
        ich_f = f["kin"]["ich"] if f.get("kin") else int(any(x != 0.0 for x in f["phirad"]))
        ich_c = c["kin"]["ich"] if c.get("kin") else int(any(x != 0.0 for x in c["phirad"]))
        if ich_f:
            n_rad += 1
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

        events.append(f["index"])
        dE_g.append(de_g)
        dp_g.append(d3_g)
        dE_e.append(de_e)
        dp_e.append(d3_e)
        dmag_e.append(dabs_e)
        rel_g.append(rg)
        rel_e.append(re_)

        mag_g = max(abs(de_g), d3_g)
        mag_e = max(abs(de_e), d3_e, abs(dabs_e))
        if mag_g > worst_g[0]:
            worst_g = (mag_g, f["index"], "photon")
        if mag_e > worst_e[0]:
            worst_e = (mag_e, f["index"], "electron")

        keep = ich_f != ich_c
        keep = keep or notable(abs(de_g), rg, args.keep_abs, args.keep_rel)
        keep = keep or notable(d3_g, d3_g / max(mag3(pg_f), mag3(pg_c), 1e-30),
                               args.keep_abs, args.keep_rel)
        keep = keep or notable(abs(de_e), re_, args.keep_abs, args.keep_rel)
        keep = keep or notable(d3_e, d3_e / max(mag3(pe_f), mag3(pe_c), 1e-30),
                               args.keep_abs, args.keep_rel)
        keep = keep or notable(abs(dabs_e), abs(dabs_e) / max(mag3(pe_f), mag3(pe_c), 1e-30),
                               args.keep_abs, args.keep_rel)
        if not keep:
            continue

        rnd = quads[f["index"] - 1] if quads[f["index"] - 1] is not None else f.get("random")
        if rnd is None:
            rnd = [float("nan")] * 4
        outliers.append({
            "event": f["index"],
            "ich_f": ich_f,
            "ich_c": ich_c,
            "r0": rnd[0],
            "r1": rnd[1],
            "r2": rnd[2],
            "r3": rnd[3],
            "vp_fx": f["vprad"][0], "vp_fy": f["vprad"][1],
            "vp_fz": f["vprad"][2], "vp_fE": f["vprad"][3],
            "vp_cx": c["vprad"][0], "vp_cy": c["vprad"][1],
            "vp_cz": c["vprad"][2], "vp_cE": c["vprad"][3],
            "pg_fx": pg_f[0], "pg_fy": pg_f[1], "pg_fz": pg_f[2], "pg_fE": pg_f[3],
            "pg_cx": pg_c[0], "pg_cy": pg_c[1], "pg_cz": pg_c[2], "pg_cE": pg_c[3],
            "pe_fx": pe_f[0], "pe_fy": pe_f[1], "pe_fz": pe_f[2], "pe_fE": pe_f[3],
            "pe_cx": pe_c[0], "pe_cy": pe_c[1], "pe_cz": pe_c[2], "pe_cE": pe_c[3],
            "dE_g": de_g,
            "dp_g": d3_g,
            "dE_e": de_e,
            "dp_e": d3_e,
            "dabs_k2": dabs_e,
            "rel_Eg": rg,
            "rel_Ee": re_,
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

    tsv = out_dir / "delta_outliers.tsv"
    fields = [
        "event", "ich_f", "ich_c", "r0", "r1", "r2", "r3",
        "vp_fx", "vp_fy", "vp_fz", "vp_fE", "vp_cx", "vp_cy", "vp_cz", "vp_cE",
        "pg_fx", "pg_fy", "pg_fz", "pg_fE", "pg_cx", "pg_cy", "pg_cz", "pg_cE",
        "pe_fx", "pe_fy", "pe_fz", "pe_fE", "pe_cx", "pe_cy", "pe_cz", "pe_cE",
        "dE_g", "dp_g", "dE_e", "dp_e", "dabs_k2", "rel_Eg", "rel_Ee",
    ]
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in outliers:
            out = {}
            for k in fields:
                v = row[k]
                out[k] = str(int(v)) if k in ("event", "ich_f", "ich_c") else fmt(v)
            w.writerow(out)

    summary_lines = [
        f"events={n}",
        f"elab={args.elab} ecm={ecm:.16e} pcm={pcm:.16e}",
        f"radiative_fortran={n_rad}",
        f"ich_mismatch={n_ich}",
        f"photon_exact_zero={n_zero_g}",
        f"electron_exact_zero={n_zero_e}",
        f"photon_nonzero={len(nonzero_g)}",
        f"electron_nonzero={len(nonzero_e)}",
        f"outliers_written={len(outliers)} keep_abs={args.keep_abs} keep_rel={args.keep_rel}",
        f"worst_photon_abs={worst_g[0]:.16e} event={worst_g[1]}",
        f"worst_electron_abs={worst_e[0]:.16e} event={worst_e[1]}",
        f"fortran={out_f}",
        f"cpp={out_c}",
        f"quads={quads_path}",
        f"outliers={tsv}",
        "photon=PHIRAD (px,py,pz,E); electron=k2 from VPRAD+PHIRAD; vp_=VPRAD",
        "deltas are C++ minus FORTRAN, GeV",
        "replay an outlier: extract_events.py quads.txt --quads --events N",
    ]
    summary = out_dir / "delta_summary.txt"
    write_summary(summary, summary_lines)
    for line in summary_lines:
        print(line)

    if not args.no_plot:
        png = out_dir / "delta_photon_electron.png"
        try:
            plot_deltas(
                png, dE_g, dp_g, dE_e, dp_e, dmag_e, events, nonzero_g, nonzero_e,
                n, n_zero_g, n_zero_e,
            )
            print(f"plot={png}")
        except ImportError:
            print("matplotlib not installed; skipped plot (data files written)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
