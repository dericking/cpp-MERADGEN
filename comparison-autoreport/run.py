#!/usr/bin/env python3
"""Run any FORTRAN vs parity-Cpp simulation and optionally write the PDF.

Dumps include HEADER (kinematics, m, m2, vpgen) and FOOTER (calls,
radiative) so the report can be rebuilt from the run directory alone.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
HARNESS = HERE / "harness"
sys.path.insert(0, str(HERE))

from lib import dump_json  # noqa: E402


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, **kw)


def write_toolchain(path: Path, build_dir: Path) -> None:
    lines = []
    for tool in ("g++", "gfortran", "cmake", "python3"):
        try:
            r = subprocess.run(
                [tool, "--version"], capture_output=True, text=True, check=True
            )
            lines.append(f"## {tool}\n{r.stdout.splitlines()[0]}\n")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            lines.append(f"## {tool}\n{e}\n")
    cache = build_dir / "CMakeCache.txt"
    if cache.is_file():
        flags = []
        for line in cache.read_text().splitlines():
            if (
                "CMAKE_BUILD_TYPE" in line
                or "CMAKE_CXX_FLAGS" in line
                or "CMAKE_Fortran_FLAGS" in line
            ):
                flags.append(line)
        lines.append("## CMakeCache (selected)\n" + "\n".join(flags) + "\n")
    path.write_text(
        "".join(lines) + "\nparity flags: -O2 -ffp-contract=off -fno-fast-math\n"
    )


def write_quads_meta(quads: Path, seed: int, n: int) -> str:
    digest = hashlib.sha256(quads.read_bytes()).hexdigest()
    quads.with_suffix(".meta").write_text(
        f"seed={seed}\nn={n}\nsha256={digest}\nfile={quads}\n"
    )
    return digest


def build(build_dir: Path) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    run(
        [
            "cmake",
            "-S",
            str(HARNESS),
            "-B",
            str(build_dir),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DPARITY_TRACE=OFF",
        ],
        check=True,
    )
    run(["cmake", "--build", str(build_dir), "-j"], check=True)


def parse_stat(path: Path) -> tuple[int, int]:
    calls = rad = -1
    for line in path.read_text().splitlines():
        if line.startswith("calls="):
            parts = line.replace("radiative=", "").split()
            calls = int(parts[0].split("=", 1)[1])
            rad = int(parts[1])
    return calls, rad


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Shared-stream FORTRAN vs meradgen-cpp run + optional PDF."
    )
    ap.add_argument("--seed", type=int, default=20260831)
    ap.add_argument("--n-radiative", type=int, default=100000)
    ap.add_argument(
        "--n-quads",
        type=int,
        default=500000,
        help="meradgen calls in the stream (must yield --n-radiative radiative events)",
    )
    ap.add_argument("--elab", type=float, default=45.0)
    ap.add_argument("--thetacm", type=float, default=90.0)
    ap.add_argument("--phi", type=float, default=10.0)
    ap.add_argument("--pl", type=float, default=-1.0)
    ap.add_argument("--prec", choices=("es14", "full", "hex"), default="full")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=HERE / "_scratch" / "run",
        help="self-contained run directory (dumps, plots, PDF)",
    )
    ap.add_argument("--skip-generate", action="store_true")
    ap.add_argument("--skip-run", action="store_true")
    ap.add_argument("--skip-analyze", action="store_true")
    ap.add_argument("--no-report", action="store_true",
                    help="stop after analyze; do not compile the PDF")
    args = ap.parse_args()

    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    build_dir = HERE / "_scratch" / "build"
    quads = (out / "quads.txt").resolve()
    out_f = (out / "fortran_radiative.txt").resolve()
    out_c = (out / "cpp_radiative.txt").resolve()

    def num(x: float) -> str:
        return str(int(x)) if float(x) == int(x) else repr(x)

    cmd = [
        "python3", "comparison-autoreport/run.py",
        "--seed", str(args.seed),
        "--n-radiative", str(args.n_radiative),
        "--n-quads", str(args.n_quads),
        "--elab", num(args.elab),
        "--thetacm", num(args.thetacm),
        "--phi", num(args.phi),
        "--pl", num(args.pl),
        "--prec", args.prec,
        "--out-dir", str(Path(args.out_dir)),
    ]

    config = {
        "date": date.today().isoformat(),
        "command": cmd,
        "seed": args.seed,
        "n_quads": args.n_quads,
        "n_radiative": args.n_radiative,
        "elab": args.elab,
        "thetacm": args.thetacm,
        "phi": args.phi,
        "pl": args.pl,
        "prec": args.prec,
        "fortran_dump": str(out_f),
        "cpp_dump": str(out_c),
        "quads": str(quads),
        "trees": {
            "fortran": "meradgen-fortran/",
            "cpp": "meradgen-cpp/",
            "cpp_final": "meradgen-cpp-final/",
        },
    }
    dump_json(out / "config.json", config)

    if not args.skip_generate:
        gen = [
            sys.executable,
            str(HARNESS / "python" / "generate_stream.py"),
            "--seed", str(args.seed),
            "--n", str(args.n_quads),
            "-o", str(quads),
        ]
        run(gen, check=True)
        digest = write_quads_meta(quads, args.seed, args.n_quads)
        print(f"quads sha256={digest}")

    if not args.skip_run:
        build(build_dir)
        fort_bin = build_dir / "meradgen_parity_fortran"
        cpp_bin = build_dir / "meradgen_parity_cpp"
        rnd = ROOT / "meradgen-fortran" / "rnd.dat"
        if rnd.is_file():
            (out / "rnd.dat").write_bytes(rnd.read_bytes())

        log_f = out / "fortran_run.log"
        log_c = out / "cpp_run.log"
        kin = [
            num(args.elab), num(args.thetacm), num(args.phi), num(args.pl),
        ]
        argv_tail = [args.prec, str(args.n_radiative), "0", *kin]
        env_f = os.environ.copy()
        env_c = os.environ.copy()
        lf = log_f.open("w")
        lc = log_c.open("w")
        try:
            pf = subprocess.Popen(
                [str(fort_bin), str(quads), str(out_f), *argv_tail],
                cwd=str(out), env=env_f, stdout=lf, stderr=subprocess.STDOUT,
            )
            pc = subprocess.Popen(
                [str(cpp_bin), str(quads), str(out_c), *argv_tail],
                env=env_c, stdout=lc, stderr=subprocess.STDOUT,
            )
            rc_f = pf.wait()
            rc_c = pc.wait()
        finally:
            lf.close()
            lc.close()
        print(f"fortran_exit={rc_f} cpp_exit={rc_c}")
        print(f"fortran_log={log_f}")
        print(log_f.read_text().strip())
        print(f"cpp_log={log_c}")
        print(log_c.read_text().strip())
        if rc_f != 0 or rc_c != 0:
            return 1
        write_toolchain(out / "toolchain.txt", build_dir)

        cf, rf = parse_stat(log_f)
        cc, rc = parse_stat(log_c)
        if rf < args.n_radiative or rc < args.n_radiative:
            print(
                f"FAIL not enough radiative fortran={rf} cpp={rc} "
                f"need={args.n_radiative} (increase --n-quads)",
                file=sys.stderr,
            )
            return 1
        if cf != cc:
            print(
                f"FAIL call-count mismatch fortran={cf} cpp={cc} "
                "(ich likely diverged)",
                file=sys.stderr,
            )
            return 1

    if not args.skip_analyze:
        rc = run([
            sys.executable, str(HERE / "analyze.py"),
            "--run-dir", str(out),
        ]).returncode
        if rc != 0:
            return rc

    if not args.no_report:
        return run([
            sys.executable, str(HERE / "write_report.py"),
            "--run-dir", str(out),
        ]).returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
