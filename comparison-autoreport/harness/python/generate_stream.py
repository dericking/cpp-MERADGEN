#!/usr/bin/env python3
"""Write N random quads using random.random() with an explicit seed (stdlib only)."""
from __future__ import annotations

import argparse
import random
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n", type=int, required=True, help="number of meradgen calls")
    ap.add_argument("--only", help="comma-separated 1-based indices to write (still draws 1..n)")
    ap.add_argument("--until", type=int, help="write events 1..until (inclusive); still draws 1..n")
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()
    if args.only and args.until is not None:
        raise SystemExit("use only one of --only / --until")
    keep = None
    if args.only:
        keep = {int(x) for x in args.only.split(",") if x.strip()}
        if not keep:
            raise SystemExit("--only needs at least one index")
        if max(keep) > args.n or min(keep) < 1:
            raise SystemExit("--only indices must be in 1..n")
    elif args.until is not None:
        if args.until < 1 or args.until > args.n:
            raise SystemExit("--until must be in 1..n")
        keep = set(range(1, args.until + 1))
    random.seed(args.seed)
    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        for i in range(1, args.n + 1):
            q = [random.random() for _ in range(4)]
            if keep is None or i in keep:
                out.write("{:.16e} {:.16e} {:.16e} {:.16e}\n".format(*q))
    finally:
        if out is not sys.stdout:
            out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
