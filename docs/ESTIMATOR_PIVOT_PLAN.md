# Estimator pivot & normalization decision — implementation plan

**Written against:** `main` @ 8 commits (post-rarefaction-removal), verified file-by-file.
**Supersedes:** `docs/PROPOSED_FIXES.md` Stages 2–4, which were written around rarefaction.
**Two independent workstreams:** Stage A (which entropy estimator) and Stage B (keep SpaNorm or not). They do not depend on each other and can run in parallel — see *Ordering* at the end.

**Status:** A0, A1, A2, B1, A3 and A4 are **done**; A5 (the estimator decision) is the next step and is deliberately left open. Sections marked ✅ record what was actually built, including the places where this plan was wrong as written — those corrections are inline below rather than in a separate errata, so the document can be followed top to bottom without cross-referencing. Everything from A3 onward is still forward-looking.

---

## Where the repo actually is

Confirmed present and working: `R/cohort.R` (D6 frozen cohort + retention floor), `R/gene_universe.R`, `R/quality_confound.R` (metric-agnostic, 418 lines), `R/entropy_correlation.R` (`pair_correlation_stats()` extracted, log10 p, Fisher-z CIs), `R/spanorm_lowmem.R`, `build_gene_universe.R`, `check_cohort_chemistry.R`, `check_cohort_retention.R`, `run_cohort.sh`, `diagnose_entropy_scaling.R`.

Confirmed removed: `sweep_rarefaction_depth.R` (folded into `rarefaction_entropy.R`), `rarefied_entropy_matrix()` / `calculate_rarefied_entropy()`, the `scuttle` dependency, the second post-exclusion depth gate, `Spots_Dropped_Entropy_Depth`.

Not yet written: `find_entropy_markers.R`, any spatial-inference module.

Three things to fix in passing, unrelated to the science:

- `environment.yml` no longer pins `scuttle`, but the archived `rarefaction_entropy.R` still calls it in three places. That script won't run in a fresh env. Either re-pin `scuttle` with a comment saying it exists only for the archived script, or add a note at the top of `rarefaction_entropy.R` saying it needs `scuttle` installed manually. Don't leave it silently broken — it's the artefact you'd reach for if an examiner asks "did you try rarefaction?"
- `README.md` still describes the pre-cohort pipeline (`analyze_entropy.R` → `analyze_stemness.R` → `find_entropy_markers.R`, an `R/` folder with two files, a `docs/` folder with one). It's badly out of date relative to its own repo. `AGENTS.md` is current. Fix README at Stage C, not before — it'll change again.
- ~~`calculate_entropy_correlations()` still defaults `entropy_cols = "shannon_entropy"`~~ ✅ **Done at A1.** The default is now `ENTROPY_COL_DEFAULT <- c("entropy_raw_plugin")`, and the `e_label` if/else chain became a named `ENTROPY_COL_LABELS` map covering the four new estimator columns, so Stage C item 4 is already half-satisfied.

Added since this plan was written:

- `tests/` — standalone check scripts run as `Rscript tests/<file>.R`, base R only, no test framework and no test-only dependency. `tests/test_entropy_estimators.R` (49 checks) is the A2 artefact.
- `R/load_sample.R` — D1's load-and-QC steps extracted out of `analyze_entropy.R` so a diagnostic can reach a QC'd Seurat object **without running SpaNorm**. See the A3 note below; this is what makes the three-sample requirement affordable.

---

## Stage A — Estimator implementation and bake-off

**Goal:** replace the two known-biased baselines with a chosen, defensible estimator, on measured evidence rather than on the argument that Chao–Shen sounds more principled.

### A0 — Freeze the success criterion *before* running anything ✅ DONE

This matters more here than it normally would, because rarefaction already taught you that a depth-correction can look successful on one number and fail on another. D6 forbids picking thresholds after seeing results; the same discipline applies to the estimator.

Frozen into `docs/DECISIONS.md` D2 before any estimator touched real data. The five criteria as written into D2 — three of them tightened from the draft below, for the reasons given:

1. **Depth decoupling.** `|ρ(H, nCount)| < 0.30` and `|ρ(H, nFeature)| < 0.30`, Spearman, on all three diagnostic samples. Not zero: depth partly tracks cellularity and more cells genuinely means more distinct programs. Set against the plug-in baseline's ρ = 0.7311 vs `nCount` and 0.8187 vs `nFeature`.
2. **No sign flip.** Concretely `ρ(H, nFeature) > −0.15` on all three samples. The draft said only "must not cross zero into moderate negative territory", which is not a threshold anyone can be held to. Rarefaction moved this statistic +0.72 → −0.23; the bound sits between those.
3. **Subsampling stability.** Median `|H(50%) − H(100%)|` under **20% of the estimator's own full-depth between-spot SD**, *and* Spearman rank correlation between the 50% and 100% values above **0.90**. The rank half is an addition: an estimator can shift level uniformly (harmless — it is monotone in the covariate of interest) or shuffle the spot ordering (fatal). The draft's single criterion cannot tell those apart. Ratio form rather than absolute bits because the estimators do not share a scale. Measured at A4.
4. **Dynamic range.** **CV = SD/mean ≥ 0.02 and IQR ≥ 0.10 bits**, per sample. The draft had no number, which makes it un-pre-registerable — "large enough to be usable" is exactly the kind of criterion that gets argued into whatever the result turned out to be. Anchored on rarefaction's documented failure (CV 0.88%) at a bit over twice it. Calibration check: the plug-in baseline on sample1 sits at CV = 2.29% (SD 0.2691, mean 11.74) and IQR = 0.2249, so the floor is live but not absurd — a candidate that fails it is genuinely worse than the metric it replaces.
5. **`percent.mt` association does not worsen.** Split in two, because the draft's version **cannot be evaluated at A3 at all**: it compares against the stemness association, and `Stemness_Score1` is produced by `analyze_stemness.R`, downstream of the estimator being chosen.
   - **5(a), at A3:** `|r(H_candidate, percent.mt)|` must not exceed `|r(H_plugin, percent.mt)|` on the same spots by more than 0.10. Self-contained inside A3's own output, and it is the form that would have caught rarefaction (r = +0.3291).
   - **5(b), at A5/Stage C:** the `percent.mt` association not larger than the entropy–stemness association, surviving the `log(nCount)` + `percent.mt` partial and holding across MT quartiles. Runs through `R/quality_confound.R` unchanged.

Criteria 1, 2, 4 and 5(a) are read off `results/statistical_tests/<sample>_estimator_comparison.csv` (A3); criterion 3 off the A4 ladder; 5(b) at Stage C.

### A1 — Implement the estimators ✅ DONE

Extend `R/shannon_entropy.R`. Everything below reuses the existing column-blocked `dgCMatrix` pattern in `plugin_entropy_matrix()` — same nnz budget, same `getOption("entropy.block.nnz", 1e6)`, same peak-memory bound. No new dependencies, consistent with D3's rejection of `ppcor` for the same reason.

**New helper: per-spot count statistics in one blocked pass.**

For each spot you need `N` (total counts), `K` (observed genes), `f1` (singletons), `f2` (doubletons). All four come off the sparse slots directly:

```r
# within a column block, after Matrix::drop0()
grp <- rep.int(seq_along(cols), diff(sub@p))
K_blk  <- diff(sub@p)
f1_blk <- tabulate(grp[sub@x == 1], nbins = length(cols))
f2_blk <- tabulate(grp[sub@x == 2], nbins = length(cols))
```

Two guards that will bite if you skip them:

- **`Matrix::drop0()` first.** Seurat's `subset()` and the gene-universe filter can leave structurally-stored explicit zeros. `diff(@p)` counts stored entries, not non-zero entries, so `K` is silently wrong without this.
- **Assert integer counts.** `f1`/`f2` are meaningless on non-integers. The raw `counts` layer should be integral; assert it rather than assume it, because this function must never be pointed at a normalized layer by accident. That mistake is the origin of the `entropy_spanorm_plugin` defect.

**Chao–Shen** (Chao & Shen 2003), the primary candidate:

- Coverage: `C = 1 - f1/N`
- Adjusted probabilities over observed genes only: `p̃_i = C · y_i / N`
- `H_CS = -Σ_i [ p̃_i · log2(p̃_i) ] / [ 1 − (1 − p̃_i)^N ]`

The denominator is a Horvitz–Thompson inclusion-probability correction: each observed gene is up-weighted by one over its probability of having been seen at all. Two numerical points:

- Compute the denominator as `-expm1(N * log1p(-p̃_i))`, not as `1 - (1-p)^N`. With `p̃ ~ 1e-4` and `N ~ 1e4` the naive form loses most of its significant digits.
- Guard `f1 == N` (every observed gene a singleton → `C = 0` → total degeneracy). Standard fix is to substitute `f1 = N - 1`. Rare at your depths but it will happen in the shallow tail, and an unguarded `NaN` propagating into a spatial plot is a bad afternoon.

**Chao–Shen with Chao–Wang–Jost coverage**, secondary candidate — same estimator, better coverage estimate using doubletons:

`C_CWJ = 1 − (f1/N) · [ (N−1)·f1 / ((N−1)·f1 + 2·f2) ]`

Worth including because `f2` is free once you're computing `f1`, and it's the natural answer if a reviewer asks "why the 2003 estimator and not something more recent." Note honestly in D2 that the *full* Chao–Wang–Jost (2013) entropy estimator is a distinct and more involved formula, and that what you implemented is Chao–Shen with the improved coverage term — don't let the label overclaim.

**Miller–Madow**, third comparator, three lines:

`H_MM = H_plugin + (K − 1) / (2 · N · ln 2)`

The `ln 2` is not optional. Your own `PROPOSED_FIXES.md` caught an earlier version of this argument stating the correction in nats while the metric was in bits, which overstated it several-fold. Getting it right here is also the cheapest possible demonstration that you understand what the correction is.

**Wrapper:** one `calculate_entropy()` Seurat wrapper taking `estimator = c("plugin", "chao_shen", "chao_shen_cwj", "miller_madow")`, writing `entropy_<estimator>`. Keep `calculate_shannon_entropy()` as-is so `diagnose_entropy_scaling.R` and the archived script don't break.

**✅ As built, with four deviations:**

- **One pass, not four.** `entropy_estimator_matrix(expr_mat, estimators = ...)` computes `N`, `K`, `f1`, `f2` and *every requested estimator* in a single blocked pass, because they all read the same four statistics. `calculate_entropy()` accepts a vector. A3 and A4 evaluate four estimators over an ~18,500 × ~11,700 matrix repeatedly; four separate passes for shared statistics is waste that compounds across the A4 ladder.
- **CWJ needs its own guards.** The plan specifies the `f1 == N` guard for Good–Turing only. The CWJ bracket `(N−1)f1 / ((N−1)f1 + 2f2)` is additionally `0/0` when `f1 = f2 = 0` — the fully-saturated spot — where the limit is a *zero* correction (`C = 1`), not `NaN`. Both guards are applied and tested.
- **Non-negativity as well as integrality.** A negative entry would pass an integrality check, produce a valid-looking `N`, and silently corrupt every proportion. Asserted together; both checks run *inside* the block loop so the temporaries stay inside the nnz budget instead of allocating a copy of the full matrix.
- **`entropy_plugin` is not aliased to `entropy_raw_plugin`.** It is bit-identical on the `counts` layer (asserted with `identical()`, not `all.equal()`), but renaming the pipeline column belongs at Stage C alongside D1's seven `*_Raw_Plugin_Entropy` QC columns. Two names on one column mid-bake-off is how a table starts lying.

Also added: `entropy_estimates(y)`, a single-count-vector convenience that routes through the production kernel rather than re-transcribing the formulas — so the A2 reference validation exercises the code the cohort run actually uses.

### A2 — Validate against a reference implementation ✅ DONE

Do not trust a hand-rolled estimator on a 11,000-spot matrix without checking it on data you can verify by hand.

Validated against CRAN **`entropy` 1.3.2** (`entropy.empirical`, `entropy.ChaoShen`, `entropy.MillerMadow`, all `unit = "log2"`) on **22 count vectors** spanning the awkward cases — uniform, all-singleton, one dominant gene plus a long tail, a single non-zero gene, a total of one count, no singletons, no doubletons, zeros interleaved, extreme skew, and Poisson / negative-binomial draws at shallow, moderate and deep sampling. **Max absolute disagreement 1.26e-12**, against the 1e-10 bar; all but one case within 2.1e-14.

Worth knowing: that residual is the *reference's*. `entropy::entropy.ChaoShen()` forms the inclusion probability as the naive `1 - (1-pa)^n` — the form this plan correctly says to avoid — so at Visium scale (`p̃ ~ 1e-4`, `N ~ 1e4`) the gap widens in the reference's disfavour. The reference is right enough to validate against at fixture scale and not the one to imitate.

Expected values frozen as literals in `tests/test_entropy_estimators.R`, `entropy` removed and confirmed absent from the pinned env. 49 checks, passing under both the conda env and system R. The file also pins the empty-spot convention (no counts → zero entropy; the reference returns `NaN`), the degeneracy guards, bitwise `entropy_plugin` ≡ `plugin_entropy_matrix()`, block-size invariance, `K`/`f1`/`f2` correctness under explicitly-stored zeros, and rejection of non-integer and negative matrices.

> **Correction — the ordering assertion in this plan was wrong.** The draft said to assert `H_plugin ≤ H_MM ≤ H_CS` spot-wise and treat any violation as a bug. Only the outer two relations are theorems. **`H_MM ≤ H_CS` is not**, and it fails on **5 of the 22 validation fixtures** — every fixture with `f1 = 0`, plus one with a light singleton tail. Miller–Madow's correction is `(K−1)/(2N ln2)`, which depends only on the observed support and keeps growing with it; Chao–Shen's is driven by `f1/N` and *vanishes* as coverage approaches 1, at which point Chao–Shen reduces to the plug-in value inflated only by the inclusion-probability denominator. So on well-sampled, low-singleton spots Miller–Madow legitimately sits above Chao–Shen.
>
> `check_entropy_ordering()` therefore **enforces** `H_plugin ≤ H_MM`, `H_plugin ≤ H_CS` and `H_plugin ≤ H_CS-CWJ` (errors on violation) and **reports** `H_MM ≤ H_CS` as an observation with a violation count. Had the full chain been asserted, A3 would have aborted on correct output. The counterexample count is itself pinned by a test, so the reasoning survives a future refactor.
>
> Related, and worth not being surprised by on the first plot: **Chao–Shen routinely exceeds `log2(K)`**, the maximum entropy of the *observed* support. That is correct — it estimates the entropy of the full distribution including unseen genes — and is asserted as expected behaviour so nobody later "fixes" it.

The real-data half — the ordering invariants spot-wise on the diagnostic samples — runs as part of A3 via `check_entropy_ordering()`.

### A3 — Depth-decoupling diagnostic ✅ DONE

Rewrite `diagnose_entropy_scaling.R` — it currently compares raw / logpac / linear-back-transform, which was the Stage 2 question and is now settled. New job: compare the estimators.

Run on **three samples spanning the depth range**: `sample1` (deep, median 17,376), `sample4` (moderate), `sample21` (shallow). Sample1 alone is what let the D = 3000 choice slip through last time; the whole point of the shallow sample is that it's where a depth-correction either works or doesn't.

Emit one tidy CSV to `results/statistical_tests/<sample>_estimator_comparison.csv` with a row per (estimator × target), targets being `nCount`, `nFeature`, `n_detected`, `percent.mt`. Reuse `pair_correlation_stats()` so the CI and log10-p conventions can't drift. Add a distribution block per estimator: mean, SD, median, IQR, range — criterion 4 is checked here.

**✅ Result summary** (full tables and interpretation in D2). Chao–Shen is unambiguously the best of the three corrections on every axis measured, and it does not clear the full A0 bar:

- **Criterion 1** — cut hard everywhere (sample21 ρ vs `nCount` 0.958 → **0.020**), but clears the 0.30 bound only on sample21; sample1 0.398 and sample4 0.578 against `nFeature` still fail.
- **Criterion 2** — passed by everything, everywhere. Nothing behaves the way rarefaction did.
- **Criterion 4** — fails on sample1 (CV 0.0126), passes on sample4 and sample21. Chao–Shen roughly halves the between-spot SD. Worth understanding rather than just recording: the correction is largest for shallow spots, so it pulls the low tail up toward the deep spots — which is simultaneously the depth decoupling it is wanted for and the range compression it is penalised for. Criteria 1 and 4 are measuring the same effect from opposite ends.
- **Criterion 5(a)** — passed everywhere, and Chao–Shen *improves* it rather than merely holding: `percent.mt` Pearson r moves −0.3701 → +0.1116 (sample1), −0.7655 → +0.1039 (sample21). Exactly the opposite of rarefaction.

Two findings that change how later stages should read: **the plug-in baseline is far worse than this project has been quoting** — ρ ≈ 0.96–0.98 vs depth on sample4 and sample21, against the sample1 figure of 0.7311 in D2 — so sample1 was flattering the baseline; and **`f1/N` is very nearly a deterministic function of depth** (ρ = −0.985 to −0.996), meaning the correction term is a monotone re-expression of depth rather than independent evidence. That is why it removes so much of the coupling and why it cannot remove all of it. **Chao–Shen-with-CWJ-coverage is not a distinguishable candidate here** — it tracks plain Chao–Shen to 3–4 decimals on every statistic; keep it as a reported comparator, don't present it as a real alternative.

**Blocker this plan didn't account for, and the fix.** Only `results/seurat_objects/sample1_spatial_obj.rds` exists. `diagnose_entropy_scaling.R` reads a saved Seurat object, so as written A3 requires `analyze_entropy.R` runs on `sample4` and `sample21` first — two SpaNorm fits at ~7.6 GB peak and hours apiece, neither of which A3 has any use for. **Every estimator in the bake-off reads raw counts only.** The single reason the old diagnostic needed a normalized object was the `logpac` layer, and that question is settled and gone from the rewrite.

So D1's load-and-QC block (steps 1–3: in-tissue subset, array coordinates, `percent.mt`/`percent.ribo`, the `nCount ≥ 500` / `nFeature ≥ 250` floors, coordinate validation, gene-universe filter) is **extracted verbatim** out of `analyze_entropy.R` into `R/load_sample.R` as `load_qc_sample()`, returning the object plus the same QC counters. `analyze_entropy.R` calls it and its output is unchanged; the diagnostic calls it and reaches an identically-QC'd object in ~1 GB and under a minute, no SpaNorm. Verified by reproducing `results/analyze_entropy/sample1_qc_metrics.csv` counter-for-counter. This also pre-pays Stage C item 1.

The diagnostic still accepts a saved `.rds` (`--from-rds`) when one exists, so sample1 can be cross-checked both ways.

### A4 — Subsampling stability (the decisive test) ✅ DONE

This is the strongest single piece of evidence you can produce, and it's the one that would have caught rarefaction's problem in advance.

Take `sample1`. Downsample every spot by a **common proportion** — 75%, 50%, 25% — and recompute all four estimators at each level. A depth-unbiased estimator returns approximately the same value at 50% as at 100%; a downward-biased one drops systematically.

Report per estimator per level: median `|ΔH|`, median signed `ΔH` (bias direction), and Spearman rank correlation against the full-depth values (does the *ranking* of spots survive, even if the level shifts).

**How the downsampling is done, now that `scuttle` is gone.** Independent binomial thinning of each stored count, `y' ~ Binomial(y, prop)`, which is dependency-free and models independent read loss. This is *not* what `scuttle::downsampleMatrix()` did — that draws without replacement per column, so each spot lands on exactly `round(prop · N)` counts, whereas thinning leaves the realized total random with mean `prop · N`. The difference is immaterial here (at `N ~ 10^4` the realized proportion varies by about 1%) and it is the right choice regardless: the point is to perturb depth, not to hit a target depth, and hitting a target depth is precisely the rarefaction assumption this project rejected. Stated so the substitution is a recorded decision rather than a silent one.

One draw per level, seeded. Averaging over draws would suppress sampling noise, but a single draw is what an experiment at that depth would actually have produced, and the *signed* median `ΔH` already separates systematic bias from noise — noise is symmetric, bias is not. That is why both the absolute and signed medians are reported rather than just the absolute one.

Two framing notes, both of which belong in D2 verbatim:

- **This uses downsampling as a measuring instrument, not as the estimator.** You are not rarefying; you are perturbing depth in a controlled way to see which estimator is invariant to it. That distinction is exactly the answer to "didn't you say rarefaction was inappropriate?" — and it's a much better answer than a citation, because it's your own data.
- Common-proportion downsampling preserves the *relative* depth structure across spots, unlike rarefaction-to-common-depth which destroys it. So this diagnostic doesn't smuggle back in the assumption you rejected.

If Chao–Shen shows near-zero median `ΔH` across the ladder while plug-in drops steadily, that is a publishable figure and the entire justification for the pivot, in one panel.

**✅ Run on all three diagnostic samples, not just sample1** — A3 had just shown sample1 to be the *least* depth-coupled of the three (plug-in ρ = 0.73 vs 0.96–0.98), so the single-sample version of this test would have been measured on the easiest case, which is the exact shape of the D = 3000 mistake this plan's closing risk warns about.

Result: **Chao–Shen removes 53–68% of the plug-in's depth bias**, and that fraction is stable across three samples and three thinning levels — a clean, reproducible finding. Miller–Madow removes 5–25%. But the figure is not the one hoped for: `ΔH` is not near zero, it is merely halved. Median signed `ΔH` at 50% depth is −0.176 / −0.163 / −0.170 bits for Chao–Shen against −0.376 / −0.490 / −0.443 for plug-in. **Criterion 3 fails for every estimator on every sample** (3–6× over the 0.20 × SD bound), and Chao–Shen's rank preservation on sample21 (ρ = 0.896) falls *below* the plug-in's (0.994) — the correction is driven by the noisy quantity `f1/N`, so it buys a smaller level shift at some cost in spot-level ranking stability. For a metric whose downstream use is a rank correlation, that trade is worth naming.

Full tables in D2. The thresholds were not revised after seeing these numbers.

### A5 — Gate and commit

Check the A0 criteria against A3/A4. Then one of:

- **Chao–Shen passes** → commit it as the primary metric. Rewrite D2 from `under review` to `implemented`, with the bake-off table as the rationale and the rejected estimators as the alternatives. This is the expected path.
- **Chao–Shen overcorrects** (sign flip, or `percent.mt` association inflates the way it did under rarefaction) → fall back to Chao–Shen-with-CWJ-coverage or Miller–Madow, which correct less aggressively. Having run all three in one pass means this costs you nothing but a paragraph.
- **Nothing decouples adequately** → the fallback is to stop trying to fix it in the *metric* and fix it in the *inference*: keep the plug-in estimator, and make the depth-adjusted partial correlation the headline number everywhere. `R/quality_confound.R` already computes exactly this. It is a weaker result but an honest one, and D3 already records the hypothesis that supports it — partialling `log(nCount)` out of the plug-in metric reproduced the rarefied metric's association (r = +0.0083 → +0.1594). If that agreement replicates, it's evidence the depth confound rather than the estimator was doing the work, and adjusting at the inference stage is legitimate.

**Effort:** A1 half a day, A2 two hours, A3 half a day including the rewrite, A4 half a day, A5 a paragraph. Call it 2 days.

---

## Stage B — Does SpaNorm earn its cost?

**Goal:** decide on evidence, not on which method sounds more sophisticated. Runs independently of Stage A; the entropy metric never touches SpaNorm either way.

Frame the question precisely, because "raw vs normalized" is the wrong frame and will lose you the argument. Entropy is a **within-spot, across-gene** statistic — only well-posed on a spot's own count proportions, which is why it goes on raw counts. The stemness score is a **between-spot, per-gene** statistic — comparing marker expression across spots, which is why it wants depth removed first. Different questions, different correct inputs, no inconsistency. Put that sentence in D1 and D4; it's the objection an examiner will raise and it's fully answerable.

What *is* open: removing a per-spot scalar size factor is one thing, and SpaNorm is doing something much stronger — a per-gene NB fit with a spatial tensor-product spline separating spatial *technical* from spatial *biological* mean trends, across ~18,500 genes. That's what costs 12.5 GB and forced `R/spanorm_lowmem.R`. The question is whether the spatial part of that is buying anything here.

### B1 — Build the hex-lattice spatial machinery ✅ DONE

Needed for Stage D anyway, so this is not extra work — just work pulled forward.

New `R/spatial_neighbors.R`, no new dependencies (this is the `spdep`-avoidance argument from `PROPOSED_FIXES.md`, and it still holds — `sf`/GDAL/GEOS/PROJ is 300–500 MB for one statistic):

- Visium spots sit on a known hex lattice. `array_col` increments by 2 within a row, so the six neighbours of `(r, c)` are `(r, c±2)`, `(r−1, c±1)`, `(r+1, c±1)`. Both `array_row` and `array_col` are already in metadata (D1) — this was foresight, use it.
- Sparse adjacency, row-standardised. ~~Under row-standardisation Moran's I reduces to `z'Wz / z'z`.~~ **Correction: only when every spot has at least one neighbour.** Under row-standardisation `S0` equals the number of *non-isolated* spots, so the general `I = (n/S0) · z'Wz / z'z` collapses to the shortcut only on an intact lattice — and the D1 QC floors punch holes in it. The factor is carried explicitly; it cancels in the permutation p-value but not in the reported I. A test pins this against a graph with an isolated spot, where the shortcut gives the wrong answer.
- Permutation p-value; 999 permutations is plenty. Two-sided `(n_extreme + 1)/(n_perm + 1)`, so it can never return 0 and never claims more evidence than the permutations support.

Built as `R/spatial_neighbors.R` and documented as **D7** (`draft` — machinery implemented and tested, no result produced with it yet). Validated by `tests/test_spatial_neighbors.R`, 33 checks, against an independent naive O(n²) transcription of the Moran's I formula sharing no code with the sparse implementation. The lattice tests assert the six offsets explicitly rather than only checking degree 6 (a transposed offset set would still give degree 6), and assert that a distant spot stays isolated rather than being attached to its nearest neighbour — the specific failure a kNN graph exhibits. One bug found and fixed: `Matrix::Diagonal() %*% adj` silently drops dimnames.

### B2 — Spatial structure in the technical covariates

Run Moran's I on `nCount_Spatial`, `nFeature_Spatial`, `percent.mt`, plus the chosen entropy metric and `Stemness_Score1`, on the three diagnostic samples.

Read it as **suggestive, not decisive** — and say so in the write-up. High spatial autocorrelation in `nCount` does *not* prove the variation is technical, because depth partly tracks cellularity and cellularity is obviously spatially structured. What it does tell you is whether there's spatially-structured variation of the kind SpaNorm's spline could act on at all. If `percent.mt` in particular shows a smooth slide-scale gradient, that's a permeabilization-style artefact and a genuine point in SpaNorm's favour.

### B3 — Does SpaNorm change the answer? (the decisive test)

Cheap, direct, and it settles the question empirically.

On the three diagnostic samples, compute the stemness module score **twice**:

1. From SpaNorm's `logpac` `data` layer — what you have now.
2. From plain Seurat `NormalizeData()` (library-size scaling + `log1p`) on the same filtered counts.

Then compare at three levels, in increasing order of what actually matters:

- **Score agreement:** Pearson and Spearman between the two `Stemness_Score1` vectors, plus side-by-side spatial plots.
- **Downstream association:** `cor(entropy_chosen, score_spanorm)` vs `cor(entropy_chosen, score_lognorm)`.
- **Conclusion stability:** run the full `run_quality_confound_check()` under both. Does the association survive the `percent.mt` + `log(nCount)` partial in both cases? Does it hold across MT quartiles in both?

The decision rule is the third one, not the first. If the *conclusion* is the same under both normalizations, SpaNorm is not earning 12.5 GB and a namespace patch for this purpose. Scores correlating at r = 0.97 while the conclusion flips would be the interesting failure mode — unlikely, but that's why you check the conclusion rather than the correlation.

One caveat to note in D4: `AddModuleScore` already partially self-corrects for depth via its control-gene-set subtraction, so some of this insensitivity is expected. That strengthens the "scalar normalization is sufficient" case rather than weakening it, but state it rather than letting a reviewer point it out.

### B4 — Branch

**Branch KEEP** — if B2 shows strong spatial structure in the technical covariates *and* B3 shows the conclusion is normalization-sensitive:

Keep everything as-is, but rewrite D1's rationale. Right now it says SCTransform and log-normalization "do not explicitly model spatial coordinates" — an argument from method description. Replace it with an argument from measurement: "we measured spatial autocorrelation in the technical covariates (Moran's I = X, p < Y) and confirmed the downstream conclusion is sensitive to the choice." That's a categorically stronger Methods paragraph, and you only get it by having run the test.

**Branch DROP** — if the conclusion is stable under plain log-normalization (the more likely outcome):

- Replace the SpaNorm call in `analyze_entropy.R` with `NormalizeData()`. Delete the ~60-line S4 `setMethod` block.
- Delete `R/spanorm_lowmem.R` and the `enable_spanorm_lowmem()` call. This is worth its own sentence: `assignInNamespace()` monkey-patching a Bioconductor package is a real maintenance and reproducibility liability, and "we removed a namespace patch because we no longer needed the package" is a good line, not an admission.
- Drop `bioconductor-spanorm` and `bioconductor-spatialexperiment` pins from `environment.yml`.
- Retire `entropy_spanorm_plugin`. It's a documented-defective reference column whose only job was to show the log-scale category error; once SpaNorm is gone it has no referent. Keep the *finding* in D2 (it's why the project pivoted) but stop computing the column. Update `run_quality_confound_check()`'s `entropy_cols` default and `calculate_entropy_correlations()`'s label map accordingly.
- Peak RSS should drop from ~7.6 GB to low single digits. At that point `run_cohort.sh` can run 2–3 samples in parallel, and the comment block explaining why it uses one process per sample needs rewriting. A 16-sample run stops being an overnight job.
- Rewrite D1. The new rejection rationale for SpaNorm must be measured, not asserted — and note explicitly that SpaNorm's *diagnosis* was valuable even though its output isn't used: it's what surfaced the log-scale degeneracy that started this whole revision.

**Middle option, and why to reject it:** keeping SpaNorm only for `find_entropy_markers.R` while scoring on log-normalized data. Tempting, but two normalizations in one pipeline is harder to defend than either alone, and you'd still pay the full memory cost. Note it as considered-and-rejected in D1 so nobody suggests it later.

**Effort:** B1 half a day, B2 two hours, B3 half a day, B4 half a day if dropping. Call it 1.5–2 days.

---

## Stage C — Consolidate and run the cohort

Do this **once**, after both A5 and B4 have landed. The whole point of doing A and B on three diagnostic samples is to avoid two 16-sample runs.

1. Rewire `analyze_entropy.R`: chosen estimator, chosen normalization, updated `entropy_cols` passed to `calculate_entropy_correlations()`.
2. Rename the seven QC distribution columns from `*_Raw_Plugin_Entropy` to match the chosen estimator. D1 lists them explicitly; keep it in sync or the QC table lies.
3. `analyze_stemness.R`: updated `entropy_cols` in both the correlation call and the confound check.
4. Update `R/entropy_correlation.R`'s `e_label` map for the new column names, and fix the `"shannon_entropy"` default noted at the top.
5. Update `DECISIONS.md` D1 (normalization), D2 (estimator, `under review` → `implemented`), D3 (status — the struck sample1 statistics get regenerated here), D4 (scoring layer + the within-spot/between-spot framing).
6. Rewrite `README.md` properly — structure, pipeline order, and the entropy/normalization description are all stale.
7. `AGENTS.md`: pipeline order line, and the SpaNorm memory note in the cohort-runs paragraph if it's gone.
8. `environment.yml`: dependency changes, plus the `scuttle` fix.
9. `./run_cohort.sh` over all 16, then `Rscript check_cohort_retention.R`.

Note that under the current `nCount ≥ 500` floor the retention gate should be a formality — near-100% everywhere, as D6 already predicts. If it isn't, something in the rewire is wrong. Free correctness check; don't skip it.

**Effort:** half a day of rewiring, plus the run.

---

## Stage D — Spatial inference (D7)

Unchanged in substance from `PROPOSED_FIXES.md`, and the reasoning still holds: ~11,700 spatially autocorrelated spots fed to `cor.test` as independent observations is why every p-value saturates. Pearson/Spearman remain fine as **effect sizes**; it's the significance test that's broken, because the effective n is far below the nominal n.

Infrastructure comes free from B1. Remaining work:

- Moran's I with permutation p on the chosen entropy metric, `nCount`, `nFeature`, `percent.mt`, `Stemness_Score1`. Partly already done in B2 — extend to the cohort.
- Spatial block permutation for the bivariate entropy–stemness question: permute in contiguous blocks so the null preserves autocorrelation.
- **Carry `percent.mt` as a covariate, not a check.** This survives the pivot: every headline number should be the partial, and the spatial model needs `percent.mt` in it. Whether the estimator is Chao–Shen or plug-in doesn't change that degradation flattens expression distributions — that's a property of the tissue, not of the estimator.

Number this **D7** in `DECISIONS.md`; D6 is taken by cohort inclusion/exclusion.

**Effort:** half a day on top of B1.

---

## Stage E — `find_entropy_markers.R`

Still unwritten, still the last step, now with one dependency worth stating: it consumes whichever normalized layer Stage B settles on. Don't start it before B4.

---

## Ordering

```
A0 (freeze criteria)  ──┐
                        ├─→ A1 → A2 → A3 → A4 → A5 ──┐
B1 (hex lattice) ───────┤                             ├─→ C → D → E
                        └─→ B2 → B3 ──────────────────┘
                                  ↑
                        (B3's "conclusion stability" half
                         needs A5's chosen metric)
```

A and B are genuinely parallel up to B3. If you want a single thread: A0 → A1 → A2 → B1 → A3 → A4 → B2 → A5 → B3 → B4 → C → D.

**Total before the cohort run: roughly 4–5 working days.**

---

## Risks worth naming upfront

**Chao–Shen may not fully decouple.** It corrects the plug-in estimator's known asymptotic downward bias; that is not the same guarantee as "no residual correlation with depth on this dataset." Your spots are heavily undersampled relative to an 18,500-gene universe, so the correction should be substantial — but substantial isn't sufficient, and A0's criteria exist precisely so you find out before committing.

**Chao–Shen could overcorrect, the same way rarefaction did.** The correction magnitude is driven by `f1/N`, which is itself depth-dependent — shallow spots get corrected harder. That's the right direction, but it's the same structural setup that produced rarefaction's sign flip. Criterion 2 exists for this. Watch it specifically.

**The `percent.mt` confound will not be fixed by any of this.** Degradation flattens the non-MT expression distribution, which inflates any diversity estimator. That's tissue biology meeting a technical artefact, and no bias-correction touches it. `R/quality_confound.R` stays a standing gate, exactly as D3 already frames it.

**One sample is not the cohort.** Every effect size in this plan gets measured on three diagnostic samples and only confirmed on sixteen. The old D = 3000 choice was made on sample1 alone and had to be revisited; don't repeat the shape of that mistake with an estimator.
