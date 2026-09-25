# Approach A — `sample_reference` / `weight_at`

MolPol needs **one** hard+soft kinematic draw (shared spectrometer tracks)
and helicity dependence only in **weights**. Today `meradgen(ppl, …)` both
samples and weights at the same `pl`, so `generate_pair(+1,−1)` can produce
different `(v,t1,z)`, `ich`, and photons.

## What depends on `pl` vs pure kinematics

| Quantity | `pl`? | Notes |
|---|---|---|
| Mandelstam `t` from `vpgen` | no | Fixed by MolPol CM angle |
| Soft+virtual `sinonr` | **yes** | `sig`, `xsBt`, `dcanc`, soft `fsirv` |
| Hard CDFs (`distsiv/t1/z`) | **yes** | built from `fsir(..., pl, ikey=2/1/0)` |
| Channel `ich` (soft vs hard) | **yes** | `P(soft)=sinonr/sitot` |
| `(v,t1,z)` draw | **yes** | inverse-CDF of those grids |
| `vectrec` → `vprad`/`phirad` | no | only `(v,t1,z)`, `vpgen`, `rand4[3]` |
| MERADGEN `weight=sitot/xs0` | **yes** | **constant in `(v,t1,z)`** at fixed `t` |

So freezing radiation is: sample once at `pl_ref`, keep 4-vectors / `(v,t1,z)` /
`ich`, and re-evaluate only polarized **densities** at that point.

## Sampling density (MERADGEN MC)

At fixed Born `t` (fixed `vpgen`):

\[
\Sigma(P)=\texttt{sitot}(P)=\texttt{sinonr}(P)+\texttt{sirad}(P),
\qquad
\texttt{sirad}(P)=\int\mathrm{d}v\,\texttt{fsir}(t,0,v,0,P;\texttt{ikey}=2).
\]

- Soft (`ich=0`): point mass at \(v=z=0\), \(t_1=t\), photon \(=0\), with
  probability \(\texttt{sinonr}/\Sigma\). Density numerator: \(\texttt{sinonr}\).
- Hard (`ich=1`): continuous density
  \[
  p_P(v,t_1,z)=\frac{1}{\Sigma(P)}\,
  \texttt{fsir}(t,t_1,v,z,P;\texttt{ikey}=0)
  \]
  (up to CDF linear-interpolation error inside grid bins).

## Unbiased importance sampling \(P_{\mathrm{ref}}\to P\)

Sample \(x\sim p_{\mathrm{ref}}\). For any test function \(f\),

\[
\int f(x)\,\sigma_P(x)\,\mathrm{d}x
=\mathbb{E}_{\mathrm{ref}}\!\left[f(x)\,\frac{\sigma_P(x)}{\sigma_{\mathrm{ref}}(x)}\right]
\cdot\Sigma_{\mathrm{ref}}.
\]

Per-event **likelihood ratio** (density numerators cancel the common \(1/\Sigma\)):

| Channel | \(\mathrm{LR}(x;P)\) |
|---|---|
| soft | \(\texttt{sinonr}(P)/\texttt{sinonr}(P_{\mathrm{ref}})\) |
| hard | \(\texttt{fsir}(\ldots P;0)/\texttt{fsir}(\ldots P_{\mathrm{ref}};0)\) |

Recommended MolPol weight (matches MERADGEN’s mean `sitot/xs0` at target \(P\)):

\[
W(x;P)=\mathrm{LR}(x;P)\cdot\frac{\Sigma(P_{\mathrm{ref}})}{\sigma_{\mathrm{Born}}(t,P)}
=\mathrm{LR}(x;P)\cdot\frac{\texttt{sitot}_{\mathrm{ref}}}{\texttt{xs0}(P)}.
\]

Identity check:

\[
\mathbb{E}_{\mathrm{ref}}[W(\cdot;P)]=\frac{\Sigma(P)}{\texttt{xs0}(P)},
\]

which is exactly `generate(P).weight` at the same `t` (constant in phase space).

Default: \(P_{\mathrm{ref}}=0\), then `weight_at(+1)` and `weight_at(-1)` on the
frozen event. Asymmetry from shared tracks:

\[
A\sim\frac{W_+-W_-}{W_++W_-}.
\]

## API (this tree)

```cpp
bool sample_reference(vp, rand4, kin, /*pl_ref=*/0);
bool weight_at(pl, vp, kin, pl_ref, WeightPieces& out);
// or low-level:
bool weight_at(pl, vp, v, t1, z, ich, sitot_ref, dens_ref, out);
```

`WeightPieces::weight` is \(W\) above. `lr` is the bare \(\mathrm{LR}\).
`sample_reference` is `generate(pl_ref, …)` — CDFs / `vectrec` run once.

## Caveats

1. **Support.** If \(\sigma_{\mathrm{ref}}(x)\approx 0\) but \(\sigma_P(x)\neq 0\),
   LR is undefined / infinite. `weight_at` returns `false`. Prefer
   \(P_{\mathrm{ref}}=0\) (symmetric support for \(\pm1\)).

2. **Soft (`ich=0`).** All soft events share identical kinematics; the LR is
   purely \(\texttt{sinonr}(P)/\texttt{sinonr}(P_{\mathrm{ref}})\). Computing
   `sinonr` needs the soft Simpson integral (same cost as MERADGEN’s soft init).
   For hard events `WeightPieces::sinonr` is left at 0 (not needed for the LR).

3. **Hard differential.** Uses `fsir(..., ikey=0)` after `zd(t,t1,v)`. Does
   **not** rebuild \(v\)/\(t_1\)/\(z\) grids or re-run `vectrec`.

4. **Grid interpolation.** The MC draws via piecewise-linear CDFs; the LR
   uses the continuous `fsir`. Residual bias is at the level of bin
   interpolation (same order as MERADGEN’s own sampling approx).

5. **Negative `fsir`.** Upstream can return negative hard densities (`nn`
   counter). LR can be negative; do not clamp. Flag / cut offline if needed.

6. **Globals.** `weight_at` mutates `t`, `pl`, `zd` coeffs, `xs0_save`, etc.
   Not thread-safe; do not interleave with a concurrent `meradgen` call.

7. **What this is not.** `WeightPieces` does not recompute integrated
   `sirad(P)` or a full `sitot(P)` per event (only via the expectation of
   `weight`). Closure tests belong in `kinematics_examination/`.

## Completeness

| Piece | Status |
|---|---|
| Sample at `pl_ref`, freeze 4-vectors | **complete** (`sample_reference`) |
| Soft LR via `sinonr` | **complete** |
| Hard LR via `fsir` ikey=0 | **complete** |
| Recommended \(W\) with \(\mathbb{E}[W]=\texttt{sitot}/\texttt{xs0}\) | **complete** |
| Integrated `sirad`/`sitot` at target `pl` without MC | not provided (not needed for event weights) |
