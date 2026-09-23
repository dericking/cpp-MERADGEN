"""Shared helpers for comparison-autoreport (stdlib only)."""
from __future__ import annotations

import csv
import json
import math
import struct
from pathlib import Path

HARNESS_PY = Path(__file__).resolve().parent / "harness" / "python"

COMPONENTS = [
    (r"VPRAD $p_x$", "vp_fx", "vp_cx"),
    (r"VPRAD $p_y$", "vp_fy", "vp_cy"),
    (r"VPRAD $p_z$", "vp_fz", "vp_cz"),
    (r"VPRAD $E$", "vp_fE", "vp_cE"),
    (r"PHIRAD $p_x$", "pg_fx", "pg_cx"),
    (r"PHIRAD $p_y$", "pg_fy", "pg_cy"),
    (r"PHIRAD $p_z$", "pg_fz", "pg_cz"),
    (r"PHIRAD $E$", "pg_fE", "pg_cE"),
]

OUTLIER_FIELDS = [
    "event", "ich_f", "ich_c", "r0", "r1", "r2", "r3",
    "vp_fx", "vp_fy", "vp_fz", "vp_fE", "vp_cx", "vp_cy", "vp_cz", "vp_cE",
    "pg_fx", "pg_fy", "pg_fz", "pg_fE", "pg_cx", "pg_cy", "pg_cz", "pg_cE",
    "pe_fx", "pe_fy", "pe_fz", "pe_fE", "pe_cx", "pe_cy", "pe_cz", "pe_cE",
    "dE_g", "dp_g", "dE_e", "dp_e", "dabs_k2", "rel_Eg", "rel_Ee",
]


def f32(x: float) -> float:
    return struct.unpack("f", struct.pack("f", float(x)))[0]


def bits32(x: float) -> int:
    return struct.unpack("I", struct.pack("f", f32(x)))[0]


def hex32(x: float) -> str:
    return f"0x{bits32(x):08X}"


def ulp32(a: float, b: float) -> int:
    if f32(a) == f32(b):
        return 0
    return abs(bits32(a) - bits32(b))


def tex_sci(x: float) -> str:
    if x == 0.0:
        return "$0$"
    mant, exp = f"{x:.3e}".split("e")
    return f"${mant}(10^{{{int(exp)}}})$"


def tex_int(n: int) -> str:
    return f"{n:,}".replace(",", "{,}")


def tex_pct(frac: float) -> str:
    return f"{100.0 * frac:.1f}\\%"


def fmt(x: float) -> str:
    return f"{x:.16e}"


def dump_diff(row: dict) -> tuple[str, float, float, int, float]:
    for label, fk, ck in COMPONENTS:
        a, b = row[fk], row[ck]
        if f32(a) != f32(b):
            return label, a, b, ulp32(a, b), abs(b - a)
    mag = max(abs(row["dE_g"]), row["dp_g"], abs(row["dE_e"]),
              row["dp_e"], abs(row["dabs_k2"]))
    return r"(reconstructed $k_2$ only)", 0.0, 0.0, 0, mag


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def dump_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def load_outliers(path: Path) -> list[dict]:
    rows = []
    if not path.is_file():
        return rows
    with path.open() as f:
        for raw in csv.DictReader(f, delimiter="\t"):
            d = {}
            for k, v in raw.items():
                if k in ("event", "ich_f", "ich_c"):
                    d[k] = int(v)
                else:
                    d[k] = float(v)
            rows.append(d)
    return rows


def write_outliers(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUTLIER_FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({
                k: (str(int(row[k])) if k in ("event", "ich_f", "ich_c")
                    else fmt(row[k]))
                for k in OUTLIER_FIELDS
            })


def load_quads(path: Path | None) -> list[list[float]]:
    if path is None or not path.is_file():
        return []
    rows = []
    with path.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                rows.append([float(parts[0]), float(parts[1]),
                             float(parts[2]), float(parts[3])])
    return rows


def k2_from(vprad: list[float], phirad: list[float],
            ecm: float, pcm: float) -> list[float]:
    return [
        -vprad[0] - phirad[0],
        -vprad[1] - phirad[1],
        pcm - vprad[2] - phirad[2],
        ecm - vprad[3] - phirad[3],
    ]


def ecm_pcm(elab: float, m: float, m2: float) -> tuple[float, float]:
    ecm = math.sqrt(2.0 * m * (elab + m)) / 2.0
    pcm = math.sqrt(ecm * ecm - m2)
    return ecm, pcm


def abs3(a: list[float], b: list[float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def mag3(v: list[float]) -> float:
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def booktabs(colspec: str, headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        r"\noindent\begin{tabularx}{\textwidth}{" + colspec + "}",
        r"\toprule",
        " & ".join(headers) + r" \\",
        r"\midrule",
    ]
    if not rows:
        ncol = len(headers)
        lines.append(rf"\multicolumn{{{ncol}}}{{c}}{{none}} \\")
    else:
        for row in rows:
            lines.append(" & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabularx}"]
    return "\n".join(lines)


def parse_meta(path: Path) -> dict:
    meta: dict = {}
    if not path.is_file():
        return meta
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if k in ("seed", "n"):
            meta[k] = int(v)
        else:
            meta[k] = v
    return meta
