#!/usr/bin/env python3
"""Write a LaTeX parity memo from a comparison-autoreport run directory."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from lib import (  # noqa: E402
    booktabs,
    dump_diff,
    hex32,
    load_json,
    load_outliers,
    tex_int,
    tex_pct,
    tex_sci,
)

MAX_ROWS = 20
MONTHS = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
]


def pretty_date(iso: str) -> str:
    y, m, d = (int(x) for x in iso.split("-"))
    return f"{d} {MONTHS[m - 1]} {y}"


def panel_rows(outliers, pred, key_fn):
    picked = [r for r in outliers if pred(r)]
    picked.sort(key=key_fn, reverse=True)
    return picked[:MAX_ROWS], len(picked)


def showing_note(n: int, n_sample: int, sort: str = r"$\max|\Delta|$") -> str:
    verb = "event" if n == 1 else "events"
    return (
        f"{n} such {verb} in {tex_int(n_sample)} "
        f"(showing {min(n, MAX_ROWS)}). Sorted by {sort}."
    )


def tex_cmd(parts: list[str]) -> str:
    if not parts:
        return r"\detokenize{comparison-autoreport/run.py}"
    toks = [r"\detokenize{" + p + "}" for p in parts]
    if len(toks) <= 2:
        return " ".join(toks)
    lines = [" ".join(toks[:2])]
    rest = toks[2:]
    for i in range(0, len(rest), 6):
        lines.append(" ".join(rest[i:i + 6]))
    return r"\\".join(lines)


def quad_preview(vals) -> str:
    if not vals:
        return r"(first line of \texttt{quads.txt})"
    return ", ".join(f"${v:.4f}$" for v in vals)


def build_tex(summary: dict, outliers: list[dict]) -> str:
    n = int(summary["n_radiative"])
    n_calls = int(summary["n_calls"])
    n_quads = int(summary.get("n_quads") or 0)
    seed = summary.get("seed")
    frac = float(summary["radiative_fraction"])
    kin = summary["kinematics"]
    elab = float(kin["elab"])
    thetacm = float(kin["thetacm"])
    phi = float(kin["phi"])
    pl = float(kin["pl"])
    n_mis = int(summary["n_dump_mismatch"])
    phys = float(summary.get("physics_rel_bar", 1e-6))
    flags = summary.get("flags", "-O2 -ffp-contract=off -fno-fast-math")
    flags_tex = flags.replace("--", r"-{}-")
    iso = str(summary.get("date", "2026-09-15"))

    p1, n1 = panel_rows(outliers, lambda r: r["dE_g"] != 0.0, lambda r: abs(r["dE_g"]))
    p3, n3 = panel_rows(outliers, lambda r: r["dE_e"] != 0.0, lambda r: abs(r["dE_e"]))
    p2, n2 = panel_rows(
        outliers,
        lambda r: r["dE_g"] != 0.0 or r["dp_g"] != 0.0,
        lambda r: max(abs(r["dE_g"]), r["dp_g"]),
    )
    p4, n4 = panel_rows(
        outliers,
        lambda r: r["dE_e"] != 0.0 or r["dp_e"] != 0.0 or r["dabs_k2"] != 0.0,
        lambda r: max(abs(r["dE_e"]), r["dp_e"], abs(r["dabs_k2"])),
    )
    six = p4 if p4 else outliers[:MAX_ROWS]

    dir_tbl = booktabs(
        r"@{}llL@{}",
        ["Directory", "In this report", "Role"],
        [
            [
                r"\texttt{meradgen-fortran/}",
                r"\textcolor{compared}{\textbf{compared}}",
                r"Upstream FORTRAN~77 reference (byte-identical to \texttt{meradgen10.tar}).",
            ],
            [
                r"\texttt{meradgen-cpp/}",
                r"\textcolor{compared}{\textbf{compared}}",
                r"Parity Cpp port: built to match FORTRAN event-by-event for the parity analysis (same constants, 32-bit \texttt{sngl} outputs).",
            ],
            [
                r"\texttt{meradgen-cpp-final/}",
                r"\textcolor{unused}{\textbf{not used}}",
                r"Production library for MolPol/Geant4: PDG constants, double precision, no \texttt{sngl}. Expected to differ from FORTRAN.",
            ],
        ],
    )

    var_tbl = booktabs(
        r"@{}llL@{}",
        [r"Code (\texttt{meradgen-fortran})", r"Math", r"Meaning in this note"],
        [
            [r"\texttt{ich}", r"---",
             r"Channel: $0$ = no real photon, $1$ = radiative."],
            [r"\texttt{vprad}", r"$p_2-p_1$",
             r"CM four-momentum of the scattered electrons, 32-bit. Dumps print this as \texttt{VPRAD}."],
            [r"\texttt{phirad}", r"$E_\gamma$",
             r"CM four-momentum of the real photon, 32-bit. $E_\gamma$ is its energy. Dumps print this as \texttt{PHIRAD}."],
            [r"---", r"$k_2$",
             r"Outgoing beam electron, rebuilt from \texttt{vprad} and \texttt{phirad}."],
            [r"\texttt{vpgen}", r"---",
             r"Input virtual-photon four-momentum ($k_1-k_2$) in the CM."],
            [r"\texttt{vgen}, \texttt{t1gen}, \texttt{zgen}", r"---",
             r"Internal double-precision kinematics. Not the official pass/fail."],
            [r"\texttt{pl}", r"$P_L$",
             rf"Longitudinal beam polarisation (here ${pl:g}$). Argument of \texttt{{meradgen}}."],
            [r"\texttt{elab}", r"$E_{\mathrm{lab}}$",
             rf"Lab beam energy (here ${elab:g}\,\mathrm{{GeV}}$)."],
            [r"\texttt{thetacm}", r"$\theta_{\mathrm{CM}}$",
             rf"CM polar angle (here ${thetacm:g}^\circ$)."],
            [r"\texttt{phi}", r"$\phi$",
             rf"CM azimuthal angle (here ${phi:g}^\circ$)."],
            [r"\texttt{urand}", r"quad",
             r"FORTRAN RNG. The shared-stream harness replaces each \texttt{urand}(\texttt{iy}) with \texttt{r1}--\texttt{r4}."],
        ],
    )

    t1 = booktabs(
        r"@{}rYYYY@{}",
        [
            r"Stream event",
            r"$E_\gamma$ FORTRAN (GeV)",
            r"$E_\gamma$ Cpp (GeV)",
            r"$\Delta E_\gamma$ (Cpp$-$FORTRAN)",
            r"$|\Delta E_\gamma|$",
        ],
        [
            [
                str(r["event"]),
                tex_sci(r["pg_fE"]),
                tex_sci(r["pg_cE"]),
                tex_sci(r["dE_g"]),
                tex_sci(abs(r["dE_g"])),
            ]
            for r in p1
        ],
    )
    t2 = booktabs(
        r"@{}rYYYYY@{}",
        [
            r"Stream event",
            r"$E_\gamma$ FORTRAN",
            r"$E_\gamma$ Cpp",
            r"$\Delta E_\gamma$",
            r"$|\Delta\vec{p}_\gamma|$",
            r"$\max|\Delta|$",
        ],
        [
            [
                str(r["event"]),
                tex_sci(r["pg_fE"]),
                tex_sci(r["pg_cE"]),
                tex_sci(r["dE_g"]),
                tex_sci(r["dp_g"]),
                tex_sci(max(abs(r["dE_g"]), r["dp_g"])),
            ]
            for r in p2
        ],
    )
    t3 = booktabs(
        r"@{}rYYYY@{}",
        [
            r"Stream event",
            r"$E_{k_2}$ FORTRAN (GeV)",
            r"$E_{k_2}$ Cpp (GeV)",
            r"$\Delta E_{k_2}$ (Cpp$-$FORTRAN)",
            r"$|\Delta E_{k_2}|$",
        ],
        [
            [
                str(r["event"]),
                tex_sci(r["pe_fE"]),
                tex_sci(r["pe_cE"]),
                tex_sci(r["dE_e"]),
                tex_sci(abs(r["dE_e"])),
            ]
            for r in p3
        ],
    )
    t4 = booktabs(
        r"@{}rYYYYYY@{}",
        [
            r"Stream event",
            r"$E_{k_2}$ FORTRAN",
            r"$E_{k_2}$ Cpp",
            r"$|\Delta E|$",
            r"$|\Delta\vec{k}_2|$",
            r"$|\Delta|k_2||$",
            r"$\max|\Delta|$",
        ],
        [
            [
                str(r["event"]),
                tex_sci(r["pe_fE"]),
                tex_sci(r["pe_cE"]),
                tex_sci(abs(r["dE_e"])),
                tex_sci(r["dp_e"]),
                tex_sci(abs(r["dabs_k2"])),
                tex_sci(max(abs(r["dE_e"]), r["dp_e"], abs(r["dabs_k2"]))),
            ]
            for r in p4
        ],
    )

    t5_rows = []
    hex_example = summary.get("hex_example")
    for r in six:
        label, a, b, ulps, ad = dump_diff(r)
        t5_rows.append([
            str(r["event"]),
            label,
            r"\texttt{" + hex32(a) + "}",
            r"\texttt{" + hex32(b) + "}",
            tex_sci(ad),
            str(ulps),
        ])
    t5 = booktabs(
        r"@{}rlYYYc@{}",
        [
            r"Stream event",
            r"Component",
            r"FORTRAN bits",
            r"Cpp bits",
            r"$|\Delta|$ (GeV)",
            r"ULPs",
        ],
        t5_rows,
    )

    if hex_example:
        hex_label = hex_example["component"]
        hex_f = hex_example["fortran_hex"]
        hex_c = hex_example["cpp_hex"]
        hex_ad = float(hex_example["abs_delta_GeV"])
        hex_u = int(hex_example["ulps"])
        hex_rel = float(hex_example.get("rel") or 0.0)
        hex_event = int(hex_example["event"])
    else:
        hex_label = hex_f = hex_c = "---"
        hex_ad = hex_rel = 0.0
        hex_u = hex_event = 0

    note1 = showing_note(n1, n, r"$|\Delta E_\gamma|$")
    note2 = showing_note(n2, n)
    note3 = showing_note(n3, n, r"$|\Delta E_{k_2}|$")
    note4 = showing_note(n4, n)

    n_tex = tex_int(n)
    calls_tex = tex_int(n_calls)
    quads_tex = tex_int(n_quads) if n_quads else r"(stream length)"
    seed_tex = str(seed) if seed is not None else r"(see \texttt{quads.meta})"
    frac_tex = tex_pct(frac)
    preview = quad_preview(summary.get("first_quad_preview") or summary.get("first_quad"))

    if summary.get("overlay_bins_identical_egamma") and summary.get("overlay_bins_identical_k2"):
        overlay_claim = (
            r"FORTRAN (solid) and Cpp (dashed) fill every bin with the same counts."
        )
        overlay_conc = (
            r"Overlay histograms of $E_\gamma$ and $|k_2|$ are bin-identical."
        )
    else:
        overlay_claim = (
            r"FORTRAN (solid) and Cpp (dashed) are overlaid; bin counts are recorded in "
            r"\texttt{bins\_egamma.tsv} and \texttt{bins\_k2.tsv}."
        )
        overlay_conc = (
            r"Overlay histograms of $E_\gamma$ and $|k_2|$ are compared bin by bin in the run directory."
        )

    if int(summary["ich_mismatch"]) == 0:
        ich_claim = r"The radiative channel \texttt{ich} never differs."
        ich_conc = r"The radiative channel never differs."
    else:
        n_ich = tex_int(int(summary["ich_mismatch"]))
        ich_claim = rf"The radiative channel \texttt{{ich}} differed in {n_ich} event(s)."
        ich_conc = ich_claim

    if int(summary["events_above_rel_1e-6"]) == 0:
        rel_claim = rf"No event exceeds relative $10^{{-6}}$ on any four-vector component."
    else:
        n_rel = tex_int(int(summary["events_above_rel_1e-6"]))
        rel_claim = (
            rf"{n_rel} event(s) exceed relative $10^{{-6}}$ on a four-vector component."
        )

    if n_mis == 0:
        mismatch_para = (
            rf"Table~\ref{{tab:alldiff}} lists every stream event whose dumped "
            r"\texttt{VPRAD} or \texttt{PHIRAD} 32-bit four-vector differs. There are "
            r"none: every dumped component matches bitwise."
        )
        mismatch_cap = r"Stream events with a nonzero dumped 32-bit four-vector. None in this sample."
        hex_para = ""
    else:
        mis_word = "event" if n_mis == 1 else "events"
        all_ulp1 = all(int(m["ulps"]) == 1 for m in summary.get("mismatches") or [])
        ulp_phrase = (
            r"Each mismatch is exactly one unit in the last place (ULP) of a "
            r"\emph{single} 32-bit component: the FORTRAN and Cpp values are "
            r"neighbouring numbers in the IEEE-754"
            if all_ulp1 else
            r"Each mismatch is measured in units in the last place (ULP) of a "
            r"32-bit component, in the IEEE-754"
        )
        factor_txt = ""
        if hex_rel > 0.0:
            factor = phys / hex_rel
            factor_txt = (
                rf", a factor of {factor:.0f} below that bar"
            )
        mismatch_para = (
            rf"Table~\ref{{tab:alldiff}} lists every stream event whose dumped "
            rf"\texttt{{VPRAD}} or \texttt{{PHIRAD}} 32-bit four-vector differs. There are "
            rf"{tex_int(n_mis)}. {ulp_phrase}\footnote{{IEEE-754 is the usual "
            r"computer standard for storing real numbers in a fixed number of bits. "
            r"A 32-bit value (also called float32 or binary32) can only represent "
            r"certain numbers exactly. The next larger or smaller allowed number is "
            r"one ULP away: the smallest step that encoding can take at that size.} "
            r"format. Occasional 1-ULP differences of this kind are expected. The two "
            r"compilers do not always add, multiply, or round in the same internal "
            r"order, even when the written formulas are the same, so they can land on "
            r"adjacent representable numbers. "
            + (
                rf"None of those {tex_int(n_mis)} reaches the relative $10^{{-6}}$ physics bar. "
                if int(summary["events_above_rel_1e-6"]) == 0 else
                rel_claim + " "
            )
            + (
                rf"The largest is stream event {tex_int(hex_event)}, "
                rf"{tex_sci(hex_ad)}\,$\mathrm{{GeV}}$ on {hex_label} (relative "
                rf"{tex_sci(hex_rel)}){factor_txt}. In hexadecimal\footnote{{A pattern such as \texttt{{{hex_f}}} is the 32 bits of the IEEE-754 number written in base~16. Values that differ by one in the last hex digit (here \texttt{{{hex_f}}} vs.\ \texttt{{{hex_c}}}) are neighbouring representable floats: one ULP.}} "
                rf"the bit pattern of that component is \texttt{{{hex_f}}} in FORTRAN and "
                rf"\texttt{{{hex_c}}} in Cpp, a {hex_u}-ULP step."
            )
        )
        mismatch_cap = (
            rf"All {tex_int(n_mis)} stream {mis_word} with a nonzero dumped 32-bit four-vector. "
            r"Sorted by $|\Delta|$."
        )
        hex_para = ""

    cmd_tex = tex_cmd(list(summary.get("command") or []))

    return rf"""\documentclass[11pt,letterpaper]{{article}}
\usepackage[T1]{{fontenc}}
\usepackage[utf8]{{inputenc}}
\usepackage{{lmodern}}
\usepackage{{microtype}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{amsmath,amssymb}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{array}}
\usepackage{{tabularx}}
\usepackage{{xcolor}}
\usepackage{{enumitem}}
\usepackage[font=small,labelfont=bf,labelsep=period]{{caption}}
\usepackage{{float}}
\usepackage{{fancyhdr}}
\usepackage[hidelinks]{{hyperref}}

\newcolumntype{{Y}}{{>{{\centering\arraybackslash}}X}}
\newcolumntype{{L}}{{>{{\raggedright\arraybackslash}}X}}

\definecolor{{compared}}{{HTML}}{{1B5E20}}
\definecolor{{unused}}{{HTML}}{{7B1E1E}}
\pagestyle{{fancy}}
\fancyhf{{}}
\renewcommand{{\headrulewidth}}{{0.4pt}}
\renewcommand{{\footrulewidth}}{{0pt}}
\lhead{{\small\textit{{MERADGEN parity note}}}}
\rhead{{\small\texttt{{meradgen-fortran}} vs.\ \texttt{{meradgen-cpp}}}}
\cfoot{{\thepage}}
\setlength{{\parskip}}{{0.55em}}
\setlength{{\parindent}}{{0pt}}
\setlist[enumerate]{{label=(\arabic*),leftmargin=1.8em,itemsep=0.25em,topsep=0.4em}}
\captionsetup[table]{{skip=4pt,position=top}}
\captionsetup[figure]{{skip=6pt}}
\setlength{{\tabcolsep}}{{4.5pt}}
\raggedbottom

\begin{{document}}
\thispagestyle{{fancy}}

\begin{{center}}
{{\LARGE\bfseries MERADGEN FORTRAN vs.\ Cpp}}\\[0.35em]
{{\large {n_tex} radiative events on a shared random stream}}\\[0.55em]
{{\normalsize {pretty_date(iso)}}}
\end{{center}}
\vspace{{0.4em}}

\section{{Introduction}}
MERADGEN is a M{{\o}}ller-scattering radiative-correction event generator.
This note compares the live FORTRAN~77 reference to the parity Cpp port
that is intended to reproduce it, including FORTRAN's float32
\texttt{{sngl}} four-vectors. ``Parity port'' here means we compare
FORTRAN and Cpp event-by-event for the parity analysis: same constants,
same 32-bit outputs. It is a port check, not a physics upgrade.
\texttt{{meradgen-cpp-final/}} is the production tree (PDG constants,
double precision throughout, no \texttt{{sngl}}) and is not used in this
comparison.

\begin{{center}}
\small
{dir_tbl}
\end{{center}}

\subsection{{How the {n_tex}-event sample was built}}
Each meradgen call draws four uniform random numbers (a ``quad'')\footnote{{MERADGEN consumes the four numbers in order: the first chooses the radiative channel (and, if radiative, samples $v$); the next two sample $t_1$ and $z$; the fourth sets a photon azimuthal sign. Example: call~1 reads the first line of \texttt{{quads.txt}}, four values in $[0,1)$ such as {preview}. The file is written by \texttt{{comparison-autoreport/run.py}} (here seed {seed_tex}, {quads_tex} quads).}}
and
returns either a non-radiative event ($\texttt{{ich}}=0$, no real photon)
or a radiative event ($\texttt{{ich}}=1$, a real photon). That choice is
part of the generator, not something we set by hand. To collect
\textbf{{{n_tex} radiative events}} we therefore had to run many more
calls and keep only those with $\texttt{{ich}}=1$.

Both directories read the \textbf{{same}} random stream (Python
\texttt{{random.random()}}, seed {seed_tex}, {quads_tex} quads available).
FORTRAN and Cpp each processed the stream in order, including
non-radiative calls, because meradgen keeps internal state from call to
call. They both reached the {n_tex}th radiative event on meradgen call
\textbf{{{calls_tex}}}. That is a radiative fraction of {frac_tex}. The plots
below contain only those {n_tex} radiative events.

Kinematics are the values written in the dump \texttt{{HEADER}}.
Defaults match \texttt{{meradgen-fortran/run.f}}:
\texttt{{elab}} $={elab:g}\,\mathrm{{GeV}}$,
\texttt{{thetacm}} $={thetacm:g}^\circ$,
\texttt{{phi}} $={phi:g}^\circ$,
\texttt{{pl}} $={pl:g}$.
Official compared outputs are the FORTRAN-style float32 four-vectors
\texttt{{VPRAD}} and \texttt{{PHIRAD}}. Deltas are Cpp minus FORTRAN, in
GeV. Internal \texttt{{KIN}} doubles (\texttt{{vgen}}, \texttt{{t1gen}},
\texttt{{zgen}}) are not the pass/fail: they can differ without moving the
float32 four-vectors. The comparison is event-by-event on those dumped
components, not a comparison of histogram bins.

\section{{Overlay spectra}}
Figure~\ref{{fig:overlay}} is the physics picture: radiated photon energy
$E_\gamma$ and outgoing CM electron momentum $|k_2|$, reconstructed from
\texttt{{VPRAD}} and \texttt{{PHIRAD}}. {overlay_claim} {ich_claim} {rel_claim}

\begin{{figure}}[H]
\centering
\includegraphics[width=\textwidth]{{overlay_radiative.png}}
\caption{{FORTRAN vs.\ Cpp overlay, {n_tex} radiative events, shared stream
through call {calls_tex}. Left: $E_\gamma$ (MeV). Right: $|k_2|$ (MeV).}}
\label{{fig:overlay}}
\end{{figure}}

\section{{Four-vector deltas (quad plot)}}
Figure~\ref{{fig:delta}} is Cpp minus FORTRAN on the same {n_tex}
radiative events. Read it left to right, top to bottom:
\begin{{enumerate}}
\item \textbf{{Photon $\Delta E_\gamma$.}} Histogram of the photon energy
difference. This panel is \emph{{only}} energy, not 3-momentum. A
difference at $10^{{-15}}\,\mathrm{{GeV}}$ sits in the $\Delta E=0$ bin at
this axis scale. All nonzero $\Delta E_\gamma$ events are listed in
Table~\ref{{tab:p1}} (up to 20).
\item \textbf{{Photon $\max|\Delta|$ vs.\ stream event.}} For each event,
$\max\bigl(|\Delta E_\gamma|,\,|\Delta\vec{{p}}_\gamma|\bigr)$. Points are
the events that are not identically zero. Table~\ref{{tab:p2}} lists
them (up to 20).
\item \textbf{{Electron $\Delta E_{{k_2}}$.}} Histogram of the
outgoing-electron energy difference. Same binning caveat as panel~(1).
All nonzero $\Delta E_{{k_2}}$ events are in Table~\ref{{tab:p3}}
(up to 20).
\item \textbf{{Electron $\max|\Delta|$ vs.\ stream event.}}
$\max\bigl(|\Delta E|,\,|\Delta\vec{{k}}_2|,\,|\Delta|k_2||\bigr)$.
Table~\ref{{tab:p4}} lists them (up to 20).
\end{{enumerate}}

\begin{{figure}}[H]
\centering
\includegraphics[width=\textwidth]{{delta_photon_electron.png}}
\caption{{Quad plot of Cpp $-$ FORTRAN, GeV. Horizontal axis on (2) and (4)
is the 1-based meradgen call index in \texttt{{quads.txt}}, not the
radiative serial number. Annotations count zeros of that panel's quantity.
$\Delta E$ is the electron energy difference (the physics quantity);
$\Delta|k_2|$ is the change in momentum length.
$|\Delta\vec{{k}}_2|$ is a port-check residual, not a quantity of physics
interest. Panel~(3) is $\Delta E$ only; panel~(4) takes the max of all three.}}
\label{{fig:delta}}
\end{{figure}}

\begin{{table}}[H]
\centering
\caption{{Panel (1) --- photon events with $\Delta E_\gamma\neq 0$.
Events whose photon 3-momentum differs but whose energy matches sit in
the $\Delta E=0$ bin of that histogram; they appear in Table~\ref{{tab:p2}}.
{note1}}}
\label{{tab:p1}}
\footnotesize
{t1}
\end{{table}}

\begin{{table}}[H]
\centering
\caption{{Panel (2) --- photon events with any nonzero $\Delta E_\gamma$ or $|\Delta\vec{{p}}_\gamma|$. {note2}}}
\label{{tab:p2}}
\footnotesize
{t2}
\end{{table}}

\begin{{table}}[H]
\centering
\caption{{Panel (3) --- electron events with $\Delta E_{{k_2}}\neq 0$.
Momentum-only differences sit in the $\Delta E=0$ bin of that histogram
and are in Table~\ref{{tab:p4}}.
{note3}}}
\label{{tab:p3}}
\footnotesize
{t3}
\end{{table}}

\begin{{table}}[H]
\centering
\caption{{Panel (4) --- electron events with any nonzero $|\Delta E|$, $|\Delta\vec{{k}}_2|$, or $|\Delta|k_2||$. {note4}}}
\label{{tab:p4}}
\footnotesize
{t4}
\end{{table}}

\clearpage
\section{{Conclusion}}
This is the equal comparison: \texttt{{meradgen-cpp/}} is built to match
\texttt{{meradgen-fortran/}}, including FORTRAN's constants and 32-bit
four-vectors. On a shared random stream it does so for {n_tex}
radiative events. {overlay_conc} {ich_conc}

{mismatch_para}
{hex_para}

\begin{{table}}[H]
\centering
\caption{{{mismatch_cap}}}
\label{{tab:alldiff}}
\footnotesize
{t5}
\end{{table}}

This does not necessarily certify the production tree
\texttt{{meradgen-cpp-final/}}. Another comparison, to ensure that
\texttt{{meradgen-cpp/}} and \texttt{{meradgen-cpp-final/}} agree within
computational expectations, will need to be performed.

\section{{Important information}}
This sample is reproduced by
\begin{{quote}}
\ttfamily\footnotesize
{cmd_tex}
\end{{quote}}
The shared stream is \texttt{{quads.txt}} in the same run directory
(seed {seed_tex}, {quads_tex} quads) if it is generated separately.
Dump \texttt{{HEADER}}/\texttt{{FOOTER}} record the kinematics actually
used and the call/radiative counts. Toolchain for this note:
\texttt{{g++}}/\texttt{{gfortran}} \texttt{{{flags_tex}}}.

\appendix
\section{{MERADGEN Variables (Code and Math)}}
\label{{app:symbols}}
\begin{{center}}
\small
{var_tbl}
\end{{center}}

\end{{document}}
"""


def compile_tex(tex_path: Path) -> None:
    subprocess.run(
        ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
         tex_path.name],
        cwd=tex_path.parent,
        check=True,
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Write the FORTRAN vs Cpp PDF from a run directory."
    )
    ap.add_argument("--run-dir", type=Path, required=True)
    args = ap.parse_args()
    run_dir = args.run_dir
    summary_path = run_dir / "summary.json"
    if not summary_path.is_file():
        print(f"FAIL missing {summary_path} (run analyze.py first)", file=sys.stderr)
        return 1
    for png in ("overlay_radiative.png", "delta_photon_electron.png"):
        if not (run_dir / png).is_file():
            print(f"FAIL missing {run_dir / png}", file=sys.stderr)
            return 1
    summary = load_json(summary_path)
    outliers = load_outliers(run_dir / "outliers.tsv")
    n = int(summary["n_radiative"])
    tex_name = f"MERADGEN_{n}_radiative_FORTRAN_vs_CPP.tex"
    pdf_name = tex_name.replace(".tex", ".pdf")
    tex_path = run_dir / tex_name
    pdf_path = run_dir / pdf_name
    tex_path.write_text(build_tex(summary, outliers))
    compile_tex(tex_path)
    copy = HERE / pdf_name
    shutil.copyfile(pdf_path, copy)
    print(f"wrote {tex_path}")
    print(f"wrote {pdf_path}")
    print(f"wrote {copy}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
