# Re-inserting a sign restriction in the exact Eo–Kim sampler

**Problem.** `eokim_exact_gibbs` (collapsed S / θ / σ²ω / (q,p) steps, mean path integrated out) is exact for the model as printed but has no sign restriction; its chains often disagree (R̂ up to 1.2, chain gap in Pr(low) 20–50 points) and some flip between labellings hundreds of times.

## Why the paper's restriction cannot go back in exactly
Under the unrestricted exact posterior the event "μ0*_t > 0 > μ1*_t at every date" has probability about 0.1% (factor data), 1% (GDP, Eo–Kim sample); 30–40% of months violate it (checked on 2,500 sweeps per chain, `/tmp/ek_pA.R` logic). The restricted posterior is therefore a drastically different, tiny-probability slice. Imposing it exactly needs the orthant probability Pr(A | y, S, θ, h) (no closed form; pseudo-marginal estimates would be hopeless at p ≈ 0.001).
The paper's own sampler imposes it by rejection (≤ 1,000 tries, then accepts anyway), which is not exact either.

## What can be imposed exactly: one linear inequality on the time-average path
Append Z_t = Σ_{s≤t} c'x_s to the Kalman state. The filtered law of Z_n at t = n is the full-sample posterior of c'Σx_t (Gaussian), so Pr(c' x̄ > 0 | y, S, θ, h) = Φ(m_Z / s_Z) in closed form (verified against 4,000 Carter–Kohn path draws: 0.9986 vs 0.9990, 0.7606 vs 0.7565).
The restricted ("joint truncation") posterior = unrestricted collapsed posterior × that probability; all MH ratios and the single-move S sweep use it (`kf_loglik_restr`, `sm_sweep_S_restr` in `R/eokim_kf.cpp`), and the path is drawn by Carter–Kohn with rejection until the restriction holds (acceptance probability = the same Φ, about 1 in practice).
`eokim_exact_gibbs(..., restrict = c(c0, c1, c2))`:  `c(1,-1,0)` mean(μ0*) > mean(μ1*) (labels ordered; no relabelling needed);  `c(0,-1,0)` mean(μ1*) < 0 (the analogue of the "low-regime mean < 0" option).

## Geweke "Getting it Right" (T = 40, M1 = 100k, M2 = 200k; prior simulator = unrestricted prior conditioned on the event by rejecting whole draws; `GEWEKE_EK_RESTRICT=diff|low Rscript geweke/run_geweke.R ekx ...`)
| sampler | seed 11 | seed 21 | seed 31 |
|---|---|---|---|
| unrestricted exact (earlier runs) | max abs Z 1.6 (this harness), 3.6 / 5.0 in earlier longer runs | | |
| restricted `low` | 1.4 | 3.5 | 2.0 |
| restricted `diff` | 2.8 | 3.6 | 4.2 (cutoff 4.39; 27 of 90 above 1.96) |
`low` is as clean as the unrestricted sampler; `diff` passes but is borderline (probably slower mixing, not a wrong target).

## Does it fix convergence? Not by itself (4 chains x 5,000 draws, `/tmp/restr_chains.R`)
| data | restriction | max gap in Pr(low) between chains | max MC s.e. | swapped-label share | AUC / F1 |
|---|---|---|---|---|---|
| factor (current) | none / diff / low | .36 / .37 / .33 | .087 / .079 / .075 | .05 / 0 / 0 | .98 / .79, .72, .74 |
| factor, COVID censored | none / diff / low | .25 / .40 / .37 | .056 / .087 / .092 | .001 / 0 / 0 | .98 / .71, .73, .71 |
| GDP, Eo–Kim sample | none / diff / low | .23 / .19 / .19 | .051 / .040 / .045 | .22 / 0 / 0 | .97 / .72, .73, .71 |
The restriction removes label switching completely (and costs about 15% more time) but the chains still disagree about which months are low (gaps of 20–40 points; ESS of s2eta 16–230 of 16,000). So the remaining problem is slow mixing of the single-move regime-path draw (a boundary between regimes cannot move several months at once), not the labelling.
**Next options:** block moves for S (shift an episode boundary / re-draw a segment with the collapsed likelihood), more or longer chains, tempering. `run_dfm_chains.R` has the model "Eo-Kim exact (a1<0)" (restrict = c(0,-1,0)) ready for the page's restriction switch; not run on the page data.

## Block moves for the regime path (`n_block`, `sm_block_S`)
MH block moves on the same collapsed target (50% boundary shifts of up to 6 months, 50% window flips of up to 12 months; both proposals symmetric, acceptance ratio = target ratio), n_block per sweep after the single-move sweep; ~18% acceptance.
**Geweke (T = 40, 10 block moves/sweep, shifts <= 3, windows <= 6, 3 seeds):** restricted `low`: max abs Z 2.4 / 1.7 / 2.9; unrestricted: 3.1 / 3.0 / 1.8 (cutoff 4.39) -> all pass.
**Real data (4 chains x 5,000, restricted `low`), max gap in Pr(low) between chains / max MC s.e. (points):**
| data | no blocks | 100 blocks | 400 blocks |
|---|---|---|---|
| factor, full sample | 33 / 7.5 | 23 / 4.7 | 26 / 5.5 |
| factor, COVID censored | 37 / 9.2 | 38 / 8.7 | 17 / 3.7 |
| GDP, Eo-Kim sample | 19 / 4.5 | 14 / 3.1 | **6 / 1.4** |
Run time per chain: factor 126 s -> 147 s -> 193 s; GDP 19 s -> 27 s -> 42 s. 400 blocks fixes the GDP sample (chains agree to within 6 points) and helps the censored factor; the full-sample factor is still 20+ points apart in some month (4 chains is a noisy estimate of the gap). The master driver uses EK_NBLOCK = 400.
