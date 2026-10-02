# Coates–Sujatha division-field computations

A bounded investigation, dated 2 October 2026, of infinite-order points on

- `E7 = [1,0,0,-1,-1]` (Cremona 294b1) over `Q(E7[7^infinity])`;
- `E5 = [1,0,0,-3,-3]` (Cremona 150a1) over `Q(E5[5^infinity])`.

These are Examples 4.7 and 4.8 of Coates–Sujatha, *Fine Selmer groups of elliptic
curves over p-adic Lie extensions*, Math. Ann. 331 (2005), 809–839.

**Neither infinite-tower question was resolved. No infinite-order point was
found.** The useful outcomes are finite-field rank bounds, a rank-zero result
for the full first 5-division field, and a way past a local bottleneck in the
septic 2-descent.

Read the [two-page handoff report](output/pdf/handoff.pdf) for the mathematical
conclusions and suggested continuation. Its [LaTeX source](report/handoff.tex)
is also included. No knowledge of this code is needed to read the report.

## Completed calculations

A row refers to the original curve base-changed to the indicated number field.

| Case | Defining polynomial | Contained in | Conclusion |
|---|---|---|---|
| 5a | x^5-x^4+4x^3+4x^2-x+13 | Q(E5[5]) | rank 0, unconditional |
| 5b | x^5+5x^3+5x-1 | Q(E5[25]) | rank 0, unconditional |
| 7a | x^7-18 | Q(E7[7]) | 2-Selmer dimension 1, rank at most 1, under GRH |
| 7b | x^7-14x^4+7x^3+28x^2+21x+8 | Q(E7[49]) | rank 0, under GRH |

The completed rank-zero descents also return trivial Sha[2], with the same
qualifications. Case 5a is isomorphic to Q(18^(1/5)). Case 7a has torsion Z/7
and root number -1. Its positive Mordell–Weil rank is **not proved**, even after
allowing GRH. Negative root number and positive Selmer rank alone do not settle
that issue.

Magma's default class-group bounds were used in fresh `-n` processes for the two
unconditional results. For the other two, the script explicitly calls
`SetClassGroupBounds("GRH")`. “GRH bound” may appear in diagnostic output even
in an unconditional run: the final proof phase, not that initial bound's label,
is what matters. The case 5b log records the ensuing finite prime checks.

## The full first 5-division field

Let L=Q(E5[5]) and F=Q(mu5). The submodule mu5 in E5[5] puts the mod-5 image
inside the affine group C5 semidirect C4, of order 20. The normal closure of
case 5a has degree 20, checked exactly, so it equals L. For K the quintic
subfield, the character identity

    Reg_G = 4 Ind_H^G(1) + Ind_N^G(1) - 4*1

with N=C5 and H=C4 gives

    rank E5(L) = 4 rank E5(K) + rank E5(F) - 4 rank E5(Q).

For completeness, on the identity the two sides have value 20; on a nonidentity
translation they have value 0; on every other element the five-point action has
one fixed point and the four-point quotient action has none, again giving 0.
Thus the character identity follows directly, and applies to E5(L) tensor Q.
The quintic rank is zero by descent, and the remaining ranks vanish by the
Selmer vanishing in Coates–Sujatha, Example 4.8. Hence **rank E5(L)=0,
unconditionally**. This does not exclude growth at higher levels.

## Reproduce

Tested with SageMath 10.7 and Magma 2.29-10. Magma requires a licence.
The Python SSH runner uses only the standard library. Point searches use the
NumPy supplied with Sage. Run from this repository's root.

```sh
sage sage/verify_fields.py
magma -n -b case_id:=5a magma/descent.m
magma -n -b case_id:=5b diagnostic:=1 magma/descent.m
magma -n -b case_id:=7b grh:=1 magma/descent.m
magma -n -b case_id:=7a grh:=1 local_hints:=1 selmer_only:=1 magma/descent.m
```

The last computation completed in 7.1 CPU seconds and printed
`TWO_SELMER_INVARIANTS [ 2 ]`. It deliberately stops before the harder cover
construction/search. Omitting `selmer_only:=1` attempts that next stage; no
successful result from that stage is claimed here.

To run over SSH with a hard time limit, one thread, low scheduling priority,
and a 4 GiB memory limit:

```sh
python3 run_remote.py --host YOUR_HOST --case 5a --seconds 600
python3 run_remote.py --host YOUR_HOST --case 7a --grh --local-hints --selmer-only --seconds 120
```

The remote host must provide `magma`, GNU `timeout`, `flock`, and `nice`.
The runner writes into `~/coates_sujatha_20261002`, with a lock preventing
concurrent computations launched by this runner. Launch one runner at a time;
it updates the remote script before acquiring the computation lock. It creates
logs in `data/`, overwriting a previous log with the same parameter name.
Direct Magma commands above do not themselves impose a wall-time cap.

For the bounded direct point search:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 sage sage/search_box.py --case 7a --basis reduced --radius 5 --seconds 300
```

The search sieves coordinate boxes modulo small primes, then makes exact
number-field square tests. Any reported point is checked against an explicit
torsion bound from good reductions. Exhausting a box proves nothing about
points outside it. Radius 5 uses substantial local memory; radius 1 is a small
check of the program. `--seconds` is an internal check during exact tests,
not a hard operating-system limit on initialisation or sieving.

## Why the local hints help

Case 7a's automatic 2-adic local-image search was the initial bottleneck.
In K=Q(a), a^7=18, the element a is a uniformizer at the prime over 2. The script
uses truncated formal-group x(t), for t=+/-a^k, to suggest local x-coordinates.
Magma verifies local solubility before using a suggestion. The factor 4 matches
the integral cubic X^3+X^2-16X-64 used internally by `TwoSelmerGroup`.
These hints are **not global points**. The completed local computation leaves
a global 2-Selmer group of order 2. Since there is no 2-torsion over K, this
gives the stated upper bound, subject to the class-group assumption.

The cache populated by `TwoSelmerGroup(E : Hints:=...)` is reused by
`MordellWeilShaInformation`. Increasing its `Effort` parameter should not be
assumed to produce a proportionately larger number-field point search. A
specialist continuation should obtain the reduced nontrivial 2-cover and search
it explicitly. No reduced covering equation is supplied as a completed result.

## Field inclusion certificates and next fields

`sage/verify_fields.py` checks polynomial identities, exact field isomorphisms,
and rational p-isogenies to 294b2 and 150a3. The latter curves' p-division fields
lie inside the original curves' p^2-division fields: the p-isogeny maps E[p^2]
onto E'[p]. Splitting groups of the four degree-p polynomials have order p(p-1).
Exact embeddings and isomorphisms are printed in `data/verify_fields.out`.

Two further candidate fields for E5, labelled 5c and 5d, are

    x^10 - 18*x^5 - 324
    x^10 - 54*x^5 + 324.

They lie in Q(E5[25]); their common discriminant is 51257812500000000.
Write alpha^5=18 and beta^5=u, where u^2-u-1=0. Then
b=beta-beta^(-1) satisfies the case 5b polynomial. For a primitive fifth root
zeta, its conjugate roots are b_j=zeta^j beta-zeta^(-j)/beta, and

    beta = (b_1-zeta^(-1)b_0)/(zeta-zeta^(-1)).

Thus beta lies in Q(E5'[5]), and alpha*beta and alpha*beta^2 define the proposed
fields in Q(E5[25]). The exact polynomial identities and irreducibility are
checked in the verification script. **No descent result is available for
5c or 5d.** Their presence as options in the scripts is preparation for future
work, not evidence of positive or zero rank.

## Logs and limitations

- `5a_tw1_grh0_e2_600s.out`: completed, rank [0,0], 25.79 CPU seconds. Its
  first version printed only the major Magma version (`VERSION 2`); the same
  installed 2.29-10 executable was used throughout.
- `5b_tw1_grh0_e2_600s_diagnostic.out`: completed, rank [0,0], 10.36 CPU seconds.
- `7b_tw1_grh1_e2_900s.out`: completed under GRH, rank [0,0], 17.28 CPU seconds.
- `7a_tw1_grh1_e2_120s_hints_selmer.out`: completed under GRH, Selmer group C2,
  7.1 CPU seconds. This is the definitive log for the septic upper bound.
- Other `7a` Magma logs: exploratory interruptions/timeouts or an SSH reset;
  none establishes a completed rank determination.
- `7b_tw1_grh0_e2_600s.out`: unconditional attempt interrupted during expensive
  certification; it establishes no unconditional rank result.
- `search_box_r5.out` and `search_box_reduced_r5.out`: completed searches in two
  bases, each with 19,487,171 coefficient vectors and 14 denominators. About
  57,000 exact tests per run, no infinite-order point. Earlier logs label the
  basis `INTEGRAL_BASIS`; the filename records which run used the reduced basis.
  The current script adds more modular sieve conditions but searches the same
  specified boxes. Small `search_smoke` logs check the final script.
- `simon_7a_600s.out`: Sage/PARI Simon descent stopped at 600 seconds; no final
  rank or point. The exploratory `sage/simon_trial.py` needs an external time
  limit if rerun.

Logs retain numerical output and failures. One server hostname is replaced by
`remote-host`; no mathematical data were changed. Failed searches are not
certificates of rank zero. The GRH-dependent calculations must not be described
as unconditional. No claim is made that the literature search establishes the
current open status of the questions in all literature.

## References

- [Coates–Sujatha (2005)](https://doi.org/10.1007/s00208-004-0609-z), Examples 4.7–4.8.
- [Rohrlich (2008)](https://doi.org/10.1090/S1056-3911-08-00462-1), §§4–5, Examples 1–2.
- [Coates–Fukaya–Kato–Sujatha (2010)](https://doi.org/10.1090/S1056-3911-09-00504-9), Corollary 6.5.

The later root-number/Selmer results do not supply the missing unconditional
Mordell–Weil point. See the report for the distinction.
