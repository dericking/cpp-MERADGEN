#!/usr/bin/env python3
"""Configure, build, run FORTRAN+C++ parity drivers, compare, record toolchain."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HARNESS = Path(__file__).resolve().parents[1]
VALROOT = Path(__file__).resolve().parents[2]


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def write_toolchain(path: Path, build_dir: Path) -> None:
    lines = []
    for tool in ("g++", "gfortran", "cmake", "python3"):
        try:
            r = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True)
            lines.append(f"## {tool}\n{r.stdout.splitlines()[0]}\n")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            lines.append(f"## {tool}\n{e}\n")
    cache = build_dir / "CMakeCache.txt"
    if cache.is_file():
        flags = []
        for line in cache.read_text().splitlines():
            if "CMAKE_BUILD_TYPE" in line or "CMAKE_CXX_FLAGS" in line or "CMAKE_Fortran_FLAGS" in line:
                flags.append(line)
        lines.append("## CMakeCache (selected)\n" + "\n".join(flags) + "\n")
    path.write_text("".join(lines) + "\nparity flags: -O2 -ffp-contract=off -fno-fast-math\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True, help="dated dir under validation_checks_new/")
    ap.add_argument("--quads", help="quads file (default: harness smoke fixture)")
    ap.add_argument("--seed", type=int, help="with --n, generate a stream into the campaign scratch")
    ap.add_argument("--n", type=int, help="number of calls when generating a stream")
    ap.add_argument("--prec", choices=("es14", "full", "hex"), default="hex",
                    help="driver print format; hex (default) is the official "
                         "compare: IEEE-754 binary32 of sngl outputs, "
                         "compare.py --bits. es14 remains a human dump.")
    ap.add_argument("--trace", action="store_true",
                    help="rebuild with MERADGEN_PARITY_TRACE and dump grid hex traces")
    args = ap.parse_args()

    campaign = VALROOT / args.campaign
    if not campaign.is_dir():
        raise SystemExit(f"no campaign directory: {campaign}")
    scratch = campaign / "_scratch"
    build = scratch / "build"
    scratch.mkdir(parents=True, exist_ok=True)
    build.mkdir(parents=True, exist_ok=True)

    if args.seed is not None:
        if args.n is None:
            raise SystemExit("--seed requires --n")
        quads = scratch / "quads.txt"
        run([sys.executable, str(HARNESS / "python" / "generate_stream.py"),
             "--seed", str(args.seed), "--n", str(args.n), "-o", str(quads)])
    elif args.quads:
        quads = Path(args.quads).resolve()
    else:
        quads = HARNESS / "fixtures" / "smoke_quads.txt"

    cmake_cmd = ["cmake", "-S", str(HARNESS), "-B", str(build),
                 "-DCMAKE_BUILD_TYPE=Release"]
    cmake_cmd.append("-DPARITY_TRACE=ON" if args.trace else "-DPARITY_TRACE=OFF")
    run(cmake_cmd)
    run(["cmake", "--build", str(build), "-j"])

    fort_bin = build / "meradgen_parity_fortran"
    cpp_bin = build / "meradgen_parity_cpp"
    out_f = scratch / "fortran_output.txt"
    out_c = scratch / "cpp_output.txt"
    rnd = ROOT / "meradgen-fortran" / "rnd.dat"
    if rnd.is_file():
        (scratch / "rnd.dat").write_bytes(rnd.read_bytes())
    env_f = os.environ.copy()
    env_c = os.environ.copy()
    trace_f = scratch / "fortran_trace.txt"
    trace_c = scratch / "cpp_trace.txt"
    if args.trace:
        env_f["MERADGEN_PARITY_TRACE"] = str(trace_f)
        env_c["MERADGEN_PARITY_TRACE"] = str(trace_c)
        for p in (trace_f, trace_c):
            if p.exists():
                p.unlink()
    run([str(fort_bin), str(quads), str(out_f), args.prec], cwd=str(scratch), env=env_f)
    run([str(cpp_bin), str(quads), str(out_c), args.prec], env=env_c)

    write_toolchain(scratch / "toolchain.txt", build)
    cmp_cmd = [sys.executable, str(HARNESS / "python" / "compare.py"), str(out_f), str(out_c)]
    if args.prec == "hex":
        cmp_cmd.append("--bits")
    cmp = subprocess.run(cmp_cmd)
    if args.trace:
        tr = subprocess.run(
            [sys.executable, str(HARNESS / "python" / "compare_trace.py"),
             str(trace_f), str(trace_c)]
        )
        if cmp.returncode == 0:
            return tr.returncode
    return cmp.returncode


if __name__ == "__main__":
    raise SystemExit(main())
