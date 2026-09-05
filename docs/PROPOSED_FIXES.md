## Stage 0 — Verify chemistry homogeneity (blocking, ~30 min)

New throwaway script `check_cohort_chemistry.R`. For each of the 22 h5 files, report to `results/cohort_qc/chemistry_check.csv`:

|column|why|
|---|---|
|`n_features`|must be identical across samples|
|`n_barcodes`|14,336 = CytAssist 11 mm capture area|
|`n_genes_all_zero`|genes with zero counts in every spot|
|`n_MT_nonzero`, `pct_MT`|MT probes present?|
|`n_RPSL_rows`, `n_RPSL_nonzero`|expect rows present, all zero, if probe-based|
|`total_counts`, `median_nCount_ontissue`|feeds the depth decision in Stage 3|

Also dump the h5 root attributes once per sample:

```r
library(hdf5r)
f <- H5File$new(p, mode = "r")
for (a in h5attr_names(f)) cat(a, ":", paste(h5attr(f, a), collapse = ", "), "\n")
f$close_all()
```

`chemistry_description` and `software_version` live there.

**Gate.** Proceed only if all 22 agree on `n_features` and show the same zero-gene pattern. sample1 reports 37,082 features with 103 `^RP[SL]` rows all at zero and MT at 3.83% — the expected explanation is a probe-based run emitting the full GRCh38 feature list, with ~19,000 genes structurally zero because no probe targets them. If `n_genes_all_zero` is around 19,000 in every sample, that hypothesis is confirmed and the cohort is homogeneous. If any sample disagrees, stop and re-plan; do not proceed to Stage 1.

**If confirmed, this collapses the original Step 0.** The probe panel is identical across all 22 samples by construction, so there is no cohort "union rule" to negotiate. The gene universe becomes: _all features with ≥1 count in ≥1 sample_ (i.e. drop the structurally-zero genes), then one detection threshold applied identically. `build_gene_universe.R` still gets written, but it is much simpler and its output is deterministic rather than a judgement call.

---

## Stage 1 — Cheap corrections (~2 h, no dependencies on anything else)

Do these first. They are unambiguous and they make everything downstream easier to read.

**Issue 5** — `analyze_stemness.R`: delete the `all(low_pcts < 2.0)` guard, which requires PROM1, POU5F1 _and_ NANOG all under 2% and therefore never fires (measured: 11.74 / 14.24 / 2.06). Replace with per-gene flags from the measured rates. Correct D4's rationale, which currently contradicts the detection table sitting next to it.

**Issue 8** — `analyze_entropy.R`: `filter.matrix = TRUE` filters the image, not the counts matrix, so the object retains all 14,336 grid barcodes. Subset explicitly on `in_tissue == 1` from `spatial/tissue_positions.csv` before depth QC, and rename `align_spots_to_coords()` to reflect what it actually does. QC denominators become `Raw_Spots_Grid`, `Spots_On_Tissue`, `Spots_Post_Depth_QC`, `Pct_OnTissue_Retained`. Final spot count is unchanged (11,709); only the accounting becomes honest. Correct D1's "aggressive cutoff removes necrotic centres" claim — it dropped 4 in-tissue spots.

**Issue 10** — `R/entropy_correlation.R`: store `Pearson_p` / `Spearman_p` as numeric doubles. Add `Pearson_log10_p` via `pt(..., log.p = TRUE)/log(10)` so it stops saturating at `0.00e+00`. Add `Pearson_CI_low/high` from `cor.test`, Spearman CI via Fisher-z, and `N`. Keep formatted strings for plot annotations only.

**Issue 3** — record measured MT/ribo content per sample from Stage 0 and correct D2's textbook "10–30% of total spot counts" claim. If probe-based is confirmed, state the _reason_ ribo is zero (excluded from the probe set), not just the number. The `exclude_pattern` becomes effectively MT-only; say so.

---

## Stage 2 — Diagnose issue 1 before rewriting anything (~1 afternoon)

The plan's diagnosis is right: `p_i = y_i / Σy_j` applied to SpaNorm's log-scale `logpac` output is a category error, and entropy degenerates to `log2(#detected genes)` (r = 0.9978). But it does **not** follow that rarefaction is the only fix, and the plan never separates the two things that are wrong.

Two distinct problems are tangled together:

1. **The log transform.** Fixable by back-transforming to linear scale.
2. **The depth-driven support ceiling.** A spot at 3,000 counts can only ever observe ~2,000 genes, so its plug-in entropy is capped near log2(2000) regardless of biology. No rescaling fixes this; the information was never captured. This is the real case for rarefaction, and it is stronger than the Miller–Madow bias argument the plan gives (that term is in nats not bits, and using K ≈ 15,000 instead of per-spot observed genes overstates it several-fold).

Run one diagnostic on sample1: take the existing `logcounts` layer, back-transform (`2^x - 1`), recompute entropy, and report `cor(entropy, nCount)`, `cor(entropy, nFeature)`, `cor(entropy, n_detected)` against the 0.680 / 0.872 / 0.998 baselines.

**Gate.** That single table tells you how much of the problem is (1) versus (2) and therefore how much machinery Stage 3 actually needs. Record it in D2 either way — it is the evidence that justifies whatever you do next.

---

## Stage 3 — Rarefied entropy (~1 day, conditional on Stage 2)

> **Superseded (2026-09-04).** This stage was implemented in full and then reverted: rarefaction is not the route this project follows. See D2 for why, and `rarefaction_entropy.R` for the code, which is kept standalone and out of the pipeline. The rest of this section is retained as the record of what was proposed. The parts that survived the reversal are the gene universe (D5) and the `entropy_raw_plugin` baseline; the parts that did not are the $D$ sweep, the `nCount \ge D` spot threshold and the `entropy_rarefied` column.

Assuming Stage 2 shows a substantial residual depth correlation after back-transformation:

`build_gene_universe.R` + `R/gene_universe.R` per Stage 0, simplified.

`sweep_rarefaction_depth.R`, one-off. D ∈ {3000, 5000, 8000}, driven by the `median_nCount_ontissue` column from Stage 0 — pick the grid from the real depth distribution rather than the plan's guessed values. Report per D: spot retention, pairwise Spearman between entropy vectors on shared spots, `cor(entropy, nCount)`, `cor(entropy, nFeature)`, `cor(entropy, n_detected)`.

`R/shannon_entropy.R`: add `calculate_rarefied_entropy()` using `scuttle::downsampleMatrix(counts, prop = depth/N, bycol = TRUE)`, averaged over `n_draws = 5`, `set.seed()` immediately before the draws rather than at file scope.

**Simplification over the original plan:** once D is committed, raise the spot QC threshold to `nCount ≥ D` instead of writing NAs for shallow spots. This removes the NA handling, removes the `listw`/data-vector length mismatch the plan then has to work around, and removes the `zero.policy` complication. Cost is a handful of spots that were going to be excluded from the entropy analysis anyway. Record the exact count in the QC table.

Columns produced:

- `entropy_rarefied` — primary
- `entropy_raw_plugin` — full-depth plug-in, demonstrates the depth bias


SpaNorm stays in the pipeline unchanged — it remains the input for DE, module scoring and visualisation. Only the entropy input changes. Say this explicitly in D1 so nobody reads the revision as "SpaNorm was a mistake."

Answers to have ready for a supervisor:

- Rarefaction discards data; alternatives are bias-corrected estimators (Miller–Madow, Chao–Shen, James–Stein) and coverage-based standardisation / Hill numbers (Chao & Jost 2012). Rarefaction is the most conservative and is standard for alpha-diversity comparison at unequal depth.
- "Rarefying is inadmissible" (McMurdie & Holmes 2014) is about differential abundance testing, not diversity estimation. It does not apply here.
- Depth partly tracks cellularity, and more cells genuinely means more distinct programs. Expect a residual depth correlation and do not treat it as failure. Revised success criterion: **substantially reduced, and no longer near-deterministic against detected-gene count.** Not "collapses to zero."

---

## Stage 4 — Spatial inference, minimum viable version (~half a day)

The problem is real: ~11,700 spatially autocorrelated spots fed to `cor.test` as independent observations is why the p-values saturate. But the original plan reaches for `spdep` + `SpatialPack`, which drags in the whole `sf`/GDAL/GEOS/PROJ stack (300–500 MB) onto a machine already peaking at 12.5 of 15 GB, to get Lee's L.

Do this instead, no new dependencies:

- Build the neighbour list directly from `array_row`/`array_col`. Visium spots sit on a known hex lattice; the six neighbour offsets are fixed. Sparse adjacency matrix, row-standardised.
- Moran's I via sparse matrix ops, with a permutation p. Run on `entropy_rarefied`, `nCount`, `nFeature`, `percent.mt`, `Stemness_Score1`. This is the evidence that motivates the whole correction and is a Methods figure in its own right.
- Spatial block permutation for the bivariate entropy-vs-stemness question: permute in contiguous blocks so the null preserves autocorrelation.

Add `spdep`/`SpatialPack` later **only if** a reviewer asks for Lee's L or a Dutilleul effective _n_ specifically, and log the dependency in D7 when you do.

**Numbering.** D6 was taken by cohort sample inclusion/exclusion (finalized). Spatial inference is **D7**.

Framing for the write-up: Pearson/Spearman are fine as _effect sizes_. What is broken is the significance test, because n ≈ 11,700 makes any p uninformative and the effective n is far smaller than the nominal n.

**Keep in view.** The effect sizes this section was originally written against were measured on the rarefied metric and no longer exist (D2). Whether there is an entropy–stemness effect worth testing properly is an open question again, to be re-answered on whatever estimator replaces the current baselines. The machinery here — Moran's I, block permutation — is still the right ceiling for a small effect and is independent of the metric, so it can be built against whatever column exists.

**Carry `percent.mt` as a covariate, not a check.** This survives the reversal on its own merits. Degradation operates at spot level inside samples that pass the sample-level D6 gate, and a degraded spot's flatter non-MT profile inflates any diversity estimator. Every headline number from here should be the partial, and the spatial model needs `percent.mt` in it.

---

## Stage 5 — Stemness panel sanity check (~2 h)

Before building more on `Stemness_Score1`:

- POU5F1 at 14.24% of spots and PROM1 at 11.74% are implausibly high for meningioma. OCT4 is notoriously prone to pseudogene multi-mapping (POU5F1B), and probe cross-hybridisation is a known issue in FFPE Visium. Check whether these are real before treating them as stemness signal.
- Add `_stemness_marker_contribution.csv`: correlation of each panel gene's expression with `Stemness_Score1`. This is what actually answers D4's stated question of whether the score reflects a multi-gene program or one dominant marker.
- Panel membership is subject to the universe rule, not force-included.

---

## Stage 6 — Docs (last, once the above is settled)

Rewrite D1–D4, add D5 (gene universe) and D6 (spatial inference) only for what was actually implemented. Record in Review notes what verification found, including the reproduction evidence and the refuted claims:

- D3's depth-decoupling claim. Nuance the original plan skipped: Spearman went _up_ after normalization (0.736 → 0.870 vs nCount, 0.822 → 0.893 vs nFeature) but Pearson went _down_ (0.680 → 0.654, 0.872 → 0.734). It is mixed, not uniformly refuted. Since the relationship is monotone-not-linear, Spearman is the right read and the claim is not supported — but state it accurately.
- D2's MT/RP abundance claim.
- D4's <2% detection claim.

Update `AGENTS.md` and `README.md` pipeline order to include stage 0 and the one-off sweep.






