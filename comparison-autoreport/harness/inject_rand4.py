#!/usr/bin/env python3
"""Copy meradgen-fortran/ into a build tree and inject the rand4 interface.

Does not edit meradgen-fortran/. The only physics-neutral change is replacing
URAND(iy) with r1,r2,r3,r4 arguments, matching meradgen-cpp.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n < 1:
        raise SystemExit(f"inject_rand4: missing {label!r} (count={n})")
    return text.replace(old, new, 1)


def insert_after(text: str, old: str, extra: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"inject_rand4: {label!r} count={n}, want 1")
    return text.replace(old, old + extra, 1)


def inject_grid_trace(text: str) -> str:
    """Harness-only Stage 3 dumps. Silent unless MERADGEN_PARITY_TRACE is set."""
    text = insert_after(
        text,
        "\t vmin=2.*Egmin*m\n",
        "      call ptrc_event()\n"
        "      call ptrc_d('t',t)\n"
        "      call ptrc_d('vmin',vmin)\n"
        "      call ptrc_d('vmax',vmax)\n",
        "vmin trace",
    )
    text = insert_after(
        text,
        "      distarv(0)=vvn\n",
        "      call ptrc_v(0,vvn,sin,distsiv(0),distarv(0))\n",
        "V0 trace",
    )
    text = insert_after(
        text,
        "      distarv(iv)=vvn\n",
        "      call ptrc_v(iv,vvn,sin,distsiv(iv),distarv(iv))\n",
        "Viv trace",
    )
    text = insert_after(
        text,
        "      sinonr=xs0+xsvr+xsB+xsf+xsadd\n",
        "      call ptrc_d('xs0',xs0)\n"
        "      call ptrc_d('xsvr',xsvr)\n"
        "      call ptrc_d('xsB',xsB)\n"
        "      call ptrc_d('xsF',xsF)\n"
        "      call ptrc_d('xsadd',xsadd)\n"
        "      call ptrc_d('sinonr',sinonr)\n",
        "sinonr trace",
    )
    text = insert_after(
        text,
        "      sitot=sirad+sinonr\n",
        "      call ptrc_d('sirad',sirad)\n"
        "      call ptrc_d('sitot',sitot)\n",
        "sitot trace",
    )
    text = insert_after(
        text,
        "      sirand=r1*sitot\n",
        "      call ptrc_d('sirand_v',sirand)\n",
        "sirand_v trace",
    )
    text = insert_after(
        text,
        "      vgen=0d0\n",
        "      call ptrc_d('vgen',vgen)\n",
        "vgen=0 trace",
    )
    text = insert_after(
        text,
        "      vgen=distarv(iv-1)+(distarv(iv)-distarv(iv-1))*\n"
        "     .(sirand-sinonr-distsiv(iv-1))/\n"
        "     .(distsiv(iv)-distsiv(iv-1))\n",
        "      call ptrc_i('v_bin',iv)\n"
        "      call ptrc_d('vgen',vgen)\n",
        "vgen interp trace",
    )
    text = insert_after(
        text,
        "      distsit1(0)=0d0\n"
        "      distart1(0)=tt1n\n",
        "      call ptrc_t(0,tt1n,sin,distsit1(0),distart1(0))\n",
        "T0 trace",
    )
    text = insert_after(
        text,
        "      distart1((i-1)*nt1+it1)=tt1n\n",
        "      call ptrc_t((i-1)*nt1+it1,tt1n,sin,\n"
        "     .distsit1((i-1)*nt1+it1),distart1((i-1)*nt1+it1))\n",
        "Tidx trace",
    )
    text = insert_after(
        text,
        "\tsirad=distsit1(4*nt1)\n",
        "      call ptrc_d('sirad_t1',sirad)\n",
        "sirad_t1 trace",
    )
    text = insert_after(
        text,
        "      sirand=r2*sirad\n",
        "      call ptrc_d('sirand_t1',sirand)\n",
        "sirand_t1 trace",
    )
    text = insert_after(
        text,
        "      t1gen=distart1(it1-1)+(distart1(it1)-distart1(it1-1))*\n"
        "     .(sirand-distsit1(it1-1))/\n"
        "     .(distsit1(it1)-distsit1(it1-1))\n",
        "      call ptrc_i('t1_bin',it1)\n"
        "      call ptrc_d('t1gen',t1gen)\n",
        "t1gen trace",
    )
    text = insert_after(
        text,
        "      distsiz(0)=0d0\n"
        "      distarz(0)=zzn\n",
        "      call ptrc_z(0,zzn,sin,distsiz(0),distarz(0))\n",
        "Z0 trace",
    )
    text = insert_after(
        text,
        "      distarz(iz)=zzn\n",
        "      call ptrc_z(iz,zzn,sin,distsiz(iz),distarz(iz))\n",
        "Ziz trace",
    )
    text = insert_after(
        text,
        "\tsirad=distsiz(nz)\n",
        "      call ptrc_d('sirad_z',sirad)\n",
        "sirad_z trace",
    )
    text = insert_after(
        text,
        "      sirand=r3*sirad\n",
        "      call ptrc_d('sirand_z',sirand)\n",
        "sirand_z trace",
    )
    text = insert_after(
        text,
        "      zgen=distarz(iz-1)+(distarz(iz)-distarz(iz-1))*\n"
        "     .(sirand-distsiz(iz-1))/\n"
        "     .(distsiz(iz)-distsiz(iz-1))\n",
        "      call ptrc_i('z_bin',iz)\n"
        "      call ptrc_d('zgen',zgen)\n",
        "zgen trace",
    )
    return text


def patch_meradgen10(src: Path) -> str:
    text = src.read_text(encoding="utf-8", errors="surrogateescape")
    text = replace_once(
        text,
        "      subroutine meradgen(ppl,vpgen)\n",
        "      subroutine meradgen(ppl,vpgen,r1,r2,r3,r4)\n",
        "meradgen signature",
    )
    text = replace_once(
        text,
        "      real*8 y,ppl,u0,u,vmax,vmin\n",
        "      real*8 y,ppl,u0,u,vmax,vmin,r1,r2,r3,r4\n",
        "meradgen locals",
    )
    text = replace_once(
        text,
        "      real*8 xsadd,sirand,urand\n",
        "      real*8 xsadd,sirand\n",
        "drop urand in meradgen",
    )
    text = replace_once(
        text,
        "      sirand=urand(iy)*sitot\n",
        "      sirand=r1*sitot\n",
        "r1*sitot",
    )
    text = replace_once(
        text,
        "      sirand=urand(iy)*sirad\n",
        "      sirand=r2*sirad\n",
        "r2*sirad",
    )
    text = replace_once(
        text,
        "      sirand=urand(iy)*sirad\n",
        "      sirand=r3*sirad\n",
        "r3*sirad",
    )
    text = replace_once(
        text,
        "\tcall vectrec(vpgen)\n",
        "\tcall vectrec(vpgen,r4)\n",
        "call vectrec",
    )
    text = replace_once(
        text,
        "      subroutine vectrec(vpgen)\n",
        "      subroutine vectrec(vpgen,r4)\n",
        "vectrec signature",
    )
    text = replace_once(
        text,
        "\treal*8 urand\n",
        "\treal*8 r4\n",
        "vectrec r4 decl",
    )
    text = replace_once(
        text,
        "      if(urand(iy).gt.0.5)sl8=-sl8\n",
        "      if(r4.gt.0.5d0)sl8=-sl8\n",
        "vectrec r4 test",
    )
    if "urand(iy)" in text:
        raise SystemExit("inject_rand4: leftover urand(iy) after patch")
    text = inject_grid_trace(text)
    return text


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, type=Path, help="meradgen-fortran directory")
    ap.add_argument("--dst", required=True, type=Path, help="build-tree copy destination")
    args = ap.parse_args()
    src = args.src.resolve()
    dst = args.dst.resolve()
    if not (src / "meradgen10.f").is_file():
        raise SystemExit(f"inject_rand4: no meradgen10.f under {src}")
    dst.mkdir(parents=True, exist_ok=True)
    inc_src = src / "include"
    inc_dst = dst / "include"
    if inc_dst.exists():
        shutil.rmtree(inc_dst)
    shutil.copytree(inc_src, inc_dst)
    shutil.copy2(src / "fsir.f", dst / "fsir.f")
    (dst / "meradgen10.f").write_text(
        patch_meradgen10(src / "meradgen10.f"),
        encoding="utf-8",
        errors="surrogateescape",
    )
    print(f"inject_rand4: wrote {dst / 'meradgen10.f'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
