
## What rarefaction did

| | ρ vs nCount | ρ vs nFeature | ρ vs percent.mt | ρ vs Stemness |
|---|---|---|---|---|
| Raw plug-in | **+0.722** | **+0.813** | −0.078 | +0.047 |
| Rarefied (D=3000) | **−0.234** | **−0.108** | **+0.298** | **+0.218** |

The depth coupling didn't collapse toward zero. It crossed zero and landed moderately negative. And a new covariate appeared: `percent.mt` at ρ = +0.298, r = +0.329, which is now the **largest** association in the table, larger than the stemness result you're actually interested in.

D3 already flags this as an open covariate, which is exactly right and is why I'm not more worried than I am. But it has to be resolved before Stage 4, not alongside it.

## The confound, and why it's the priority

I think the negative depth correlation and the positive MT correlation are the same artifact. The likely mechanism: a high-MT spot is degraded, so its non-MT signal is flatter and closer to ambient noise, which at fixed depth reads as *higher* diversity. And low-depth spots tend to be high-MT (cohort-level ρ = −0.41 from your own Stage 0 table). So both columns are plausibly measuring the same thing: **low-quality spots score high on rarefied entropy.**

That matters because your stemness result moved from r = +0.008 (null) to r = +0.184 in the same step that introduced the MT association. If `Stemness_Score1` also tracks MT — plausible, since `AddModuleScore` on a 7-gene panel behaves oddly in sparse degraded spots — then the entropy-stemness correlation could be entirely quality-driven.

This is roughly thirty lines of work and it decides whether Stage 4 is worth running:

1. Report `cor(Stemness_Score1, percent.mt)` and `cor(Stemness_Score1, nCount)`. These aren't in any output file yet.
2. Partial correlation of entropy against stemness, controlling for `percent.mt` and `log(nCount)`. Either `ppcor::pcor.test`, or residualize both on the covariates and correlate the residuals.
3. Subset check: restrict to spots with `percent.mt < 10` and rerun the raw correlation. If the association holds in the clean subset, it's real.

If the partial correlation survives, you have a genuine finding and Stage 4 gives it rigorous inference. If it collapses, Stage 4 would be attaching careful spatial p-values to an artifact.

## Two remaining correctness gaps

**D = 3000 hasn't been validated cohort-wide.** The sweep ran on sample1 only, which is one of your deepest samples (median 17,376). My D = 2000 recommendation was driven by sample21 and sample4, not sample1. Worse, the sweep correctly measures retention on *post-exclusion* depth, which is stricter than the `nCount` quantiles I estimated from — for your 16 retained samples at 3.4–8.7% MT, post-exclusion depth runs 91–97% of `nCount`. So sample21, which I estimated at ~71% retention at D = 3000, will come in lower.

Run `sweep_rarefaction_depth.R sample21` and `sample4` before committing. It's two commands.

The good news from the pairwise matrix: ρ between D = 2000 and D = 3000 is 0.966, and between 3000 and 8000 it's 0.982. Spot rankings are nearly invariant to D above 1000, so dropping to D = 2000 for retention costs you very little in rank-based analysis. D = 1000 is the outlier (ρ ≈ 0.94 against everything else) and should stay off the table.

One thing to watch: `SD_Rarefied_Entropy` is 0.098 against a mean of 10.61, a coefficient of variation under 1%. The metric discriminates, but the dynamic range is narrow, and it does grow with D (0.089 at 2000, 0.116 at 8000). Report the IQR and full range alongside the SD so the narrowness is visible rather than buried.

**The per-sample gene filter is still there and still breaks comparability.** Both `analyze_entropy.R` and `sweep_rarefaction_depth.R` apply `min_spots_per_gene <- max(20, ceiling(0.02 * n_spots))` *after* the universe filter. So sample1 ends at 15,256 genes and other samples will end elsewhere. A gene kept in sample2 but filtered in sample1 contributes to sample2's entropy and to its depth accounting, and doesn't in sample1. The universe fixed the starting point but not the finishing one.

Decide before running samples 2–22. Either apply one cohort-wide detection rule producing a single frozen gene set, or drop the per-sample filter entirely and use the universe as-is — rarefaction already handles the depth problem the filter was partly guarding against, and rare genes contribute little to entropy. I'd lean toward dropping it, since it's simpler to defend, but either is fine as long as it's identical across samples.

## Still missing from the docs

The 6-sample exclusion is still not written down. I flagged it last round; it's still absent. Everything so far is sample1, which survives the cut, so it hasn't bitten yet. But the moment you run samples 2–22 and see their entropy-stemness numbers, a decision recorded afterward is worth much less than one recorded before. Write it now.

The Space Ranger version discrepancy (files say 2.0.1, GEO says 4.0.1) is also still unrecorded.

## How to proceed

Stage 4 as specified is still Moran's I plus a spatial block permutation, built from the `array_row`/`array_col` hex lattice with no new dependencies. But four things come first:

1. The MT/quality confounding check on sample1. This is the gate — if the stemness association doesn't survive it, everything downstream changes.
2. Sweep sample21 and sample4, then commit D.
3. Freeze the gene filter cohort-wide.
4. Write the exclusion decision and the version note into DECISIONS.md.
5. Run all 16 samples.

Then Stage 4. Doing the spatial inference before step 1 would mean building careful machinery to test an effect whose provenance you don't yet know.
