#!/usr/bin/env python3
"""Extract selected EVENT blocks or quad rows by 1-based index. Stdlib only."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import EVENT_RE  # noqa: E402


def parse_indices(s: str) -> list[int]:
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(int(part))
    if not out:
        raise SystemExit("no indices")
    return out


def extract_output(path: Path, want: set[int]) -> str:
    lines = path.read_text().splitlines(True)
    chunks: dict[int, list[str]] = {}
    cur = None
    buf: list[str] = []
    for line in lines:
        m = EVENT_RE.match(line)
        if m:
            if cur is not None and cur in want:
                chunks[cur] = buf
            cur = int(m.group(1))
            buf = [line]
            continue
        if cur is not None:
            buf.append(line)
    if cur is not None and cur in want:
        chunks[cur] = buf
    missing = sorted(want - set(chunks))
    if missing:
        raise SystemExit(f"missing events: {missing}")
    return "".join("".join(chunks[i]) for i in sorted(want))


def extract_quads(path: Path, want: list[int]) -> str:
    rows = path.read_text().splitlines()
    out = []
    for i in want:
        if i < 1 or i > len(rows):
            raise SystemExit(f"quad index {i} out of range 1..{len(rows)}")
        out.append(rows[i - 1] + "\n")
    return "".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="event dump or quads.txt")
    ap.add_argument("--events", required=True, help="comma-separated 1-based indices")
    ap.add_argument("--quads", action="store_true", help="treat file as one quad per line")
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()
    want_list = parse_indices(args.events)
    path = Path(args.path)
    text = extract_quads(path, want_list) if args.quads else extract_output(path, set(want_list))
    if args.output == "-":
        sys.stdout.write(text)
    else:
        Path(args.output).write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
