# Independent Audit Report: Repository State & Estimator Pivot Plan

**Target Repository:** `meningioma-spatial-entropy` (`/home/alexej/Projects/Tirocinio`)  
**Audit Reference:** [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md), [DECISIONS.md](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md), [AGENTS.md](file:///home/alexej/Projects/Tirocinio/AGENTS.md)  
**Git HEAD at Review:** `4b58e5d515fbdeaef90b04a1a1101e38da091d98`  
**Date:** September 7, 2026  

---

## 1. Executive Summary

This report provides an objective, rigorous review of the repository's current state relative to the plan set forth in [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md).

### Overall Assessment
1. **Mathematical & Statistical Rigor:** The core implementations of the entropy estimators, reference validations, subsampling stability ladder, and hex-lattice spatial autocorrelation machinery are **exceptionally sound**. There are no hallucinations in the mathematical derivations (Chao–Shen Horvitz–Thompson inclusion probability, Chao–Wang–Jost coverage, Miller–Madow bits conversion, Moran's $I$ degree adjustment). Unit tests pass 100% (82/82 checks across two standalone suites).
2. **Critical Git Tracking Discrepancy:** The last two git commits (`e018f2d` and `4b58e5d`) committed **only markdown documentation** ([DECISIONS.md](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md) and [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md)). **None of the R implementation files, test scripts, diagnostic tools, or sweep drivers were committed.** The entire codebase implementation currently resides uncommitted in the working tree (as untracked files or unstaged modifications). If cloned from remote, the repository is missing all newly developed code.
3. **Progress Milestone:** The project has completed **all of Stage A (A0 through A5b)** and **Stage B1**. Stages **B2, B3, and B4** (evaluating whether SpaNorm earns its memory/time overhead vs. plain log-normalization) have **not been executed**.
4. **Premature Execution / Sequencing Fault:** The developer prematurely jumped into **Stage C** by modifying [analyze_entropy.R](file:///home/alexej/Projects/Tirocinio/analyze_entropy.R) and [analyze_stemness.R](file:///home/alexej/Projects/Tirocinio/analyze_stemness.R) and running `sample1` through SpaNorm before resolving Stage B. If Stage B decides to drop SpaNorm (which the plan anticipates as the most likely outcome), these modifications and the expensive SpaNorm execution on `sample1` will have to be redone.
5. **Minor Code / Documentation Bugs:** A minor statistical implementation flaw was identified in [diagnose_correction_routes.R](file:///home/alexej/Projects/Tirocinio/diagnose_correction_routes.R) (semi-rank correlation instead of full Spearman correlation), an internal status inconsistency exists in [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md), and a factual inconsistency was found in [DECISIONS.md](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md) regarding Criterion 4 cohort pass rates.

---

## 2. Git Commit Audit vs. Working Tree State

### Recent Commit History
```
4b58e5d (HEAD -> main, origin/main) further implementation of the new plan
e018f2d implementation of the new plan up to A4
9fe8add Merge pull request #1 from arex23/remove-rarefaction
```

### Commit Content vs. Commit Message Reality
- **Commit `e018f2d`** ("implementation of the new plan up to A4"):
  - **Files changed:** [DECISIONS.md](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md) (+117, -1), [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md) (+319, -1).
  - **Actual code committed:** **0 lines**.
- **Commit `4b58e5d`** ("further implementation of the new plan"):
  - **Files changed:** [DECISIONS.md](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md) (+152, -15), [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md) (+36, -2).
  - **Actual code committed:** **0 lines**.

### Uncommitted Working Tree Inventory
The actual implementations are stranded in the local working tree:
- **Untracked Core Files:**
  - [R/load_sample.R](file:///home/alexej/Projects/Tirocinio/R/load_sample.R) (Extracts D1 load & QC to allow fast diagnostic loading without SpaNorm)
  - [R/spatial_neighbors.R](file:///home/alexej/Projects/Tirocinio/R/spatial_neighbors.R) (Stage B1 hex-lattice graph and Moran's $I$)
  - [R/diagnostic_sweep.R](file:///home/alexej/Projects/Tirocinio/R/diagnostic_sweep.R) (Stage A cohort sweep driver)
  - [diagnose_subsampling_stability.R](file:///home/alexej/Projects/Tirocinio/diagnose_subsampling_stability.R) (Stage A4 subsampling ladder)
  - [diagnose_correction_routes.R](file:///home/alexej/Projects/Tirocinio/diagnose_correction_routes.R) (Stage A5a route comparison)
  - [test_entropy_estimators.R](file:///home/alexej/Projects/Tirocinio/tests/test_entropy_estimators.R) (Stage A2 unit tests, 49 tests)
  - [test_spatial_neighbors.R](file:///home/alexej/Projects/Tirocinio/tests/test_spatial_neighbors.R) (Stage B1 unit tests, 33 tests)
  - All cohort sweep output CSVs and PNGs under `results/statistical_tests/cohort_*`.
- **Modified Tracked Files (Unstaged):**
  - [R/shannon_entropy.R](file:///home/alexej/Projects/Tirocinio/R/shannon_entropy.R) (+355 lines: estimator kernels, CWJ coverage, Miller–Madow, `calculate_entropy()`)
  - [R/entropy_correlation.R](file:///home/alexej/Projects/Tirocinio/R/entropy_correlation.R) (+53 lines: label mapping, default columns)
  - [R/quality_confound.R](file:///home/alexej/Projects/Tirocinio/R/quality_confound.R) (+12 lines: label delegation)
  - [diagnose_entropy_scaling.R](file:///home/alexej/Projects/Tirocinio/diagnose_entropy_scaling.R) (+423, -200 lines: rewritten for A3 bake-off)
  - [AGENTS.md](file:///home/alexej/Projects/Tirocinio/AGENTS.md) (+27 lines: documentation of diagnostics)
  - [analyze_entropy.R](file:///home/alexej/Projects/Tirocinio/analyze_entropy.R) & [analyze_stemness.R](file:///home/alexej/Projects/Tirocinio/analyze_stemness.R) (Premature Stage C edits)
  - Various `results/` artifacts from running `sample1`.

> [!WARNING]
> **Reproducibility Risk:** Any external collaborator pulling `origin/main` receives documentation claiming the estimator pivot is complete and validated, but none of the scripts or functions will exist in their clone.

---

## 3. Plan Progression Milestone Analysis

| Stage | Plan Description | Code State | Data / Execution State | Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **A0** | Freeze acceptance criteria before testing | Recorded in `DECISIONS.md` D2 | Criteria frozen prior to runs | **COMPLETE** |
| **A1** | Implement entropy estimators | Implemented in [R/shannon_entropy.R](file:///home/alexej/Projects/Tirocinio/R/shannon_entropy.R) | Single-pass blocked kernel | **COMPLETE** |
| **A2** | Reference validation against CRAN `entropy` | Implemented in [test_entropy_estimators.R](file:///home/alexej/Projects/Tirocinio/tests/test_entropy_estimators.R) | 49/49 tests pass ($\Delta \le 1.26\times 10^{-12}$) | **COMPLETE** |
| **A3** | Depth-decoupling diagnostic | Implemented in [diagnose_entropy_scaling.R](file:///home/alexej/Projects/Tirocinio/diagnose_entropy_scaling.R) | Run on 3 diagnostic + 16 cohort samples | **COMPLETE** |
| **A4** | Subsampling stability ladder | Implemented in [diagnose_subsampling_stability.R](file:///home/alexej/Projects/Tirocinio/diagnose_subsampling_stability.R) | Run on 3 diagnostic + 16 cohort samples | **COMPLETE** |
| **A5** | Estimator decision gate (A5a & A5b) | Implemented in [diagnose_correction_routes.R](file:///home/alexej/Projects/Tirocinio/diagnose_correction_routes.R) | 16-sample sweep executed, decision documented | **COMPLETE** |
| **B1** | Hex-lattice spatial neighborhood & Moran's $I$ | Implemented in [R/spatial_neighbors.R](file:///home/alexej/Projects/Tirocinio/R/spatial_neighbors.R) | 33/33 tests pass in [test_spatial_neighbors.R](file:///home/alexej/Projects/Tirocinio/tests/test_spatial_neighbors.R) | **COMPLETE** |
| **B2** | Spatial autocorrelation of technical covariates | Not written | Not executed | **NOT STARTED** |
| **B3** | Decisive test: SpaNorm vs. `NormalizeData()` | Not written | Not executed | **NOT STARTED** |
| **B4** | Branch decision: KEEP or DROP SpaNorm | Undecided | Depends on B2/B3 | **NOT STARTED** |
| **C** | Pipeline consolidation & full cohort run | Partially modified `analyze_*.R` | Only `sample1` run; cohort sweep pending | **PREMATURE / PARTIAL** |
| **D** | Spatial inference (block permutations) | Not written | Not executed | **NOT STARTED** |
| **E** | Marker discovery (`find_entropy_markers.R`) | Not written | Not executed | **NOT STARTED** |

---

## 4. Evaluation of Hallucinations, Wrong Choices, and Discrepancies

### A. Code & Methodological Soundness (No Scientific Hallucinations)
A deep code review of the mathematical modules confirmed high rigor:
1. **Numerical Stability:** The Horvitz–Thompson inclusion probability is implemented as:
   $$\lambda = -\text{expm1}(N \cdot \text{log1p}(-\tilde{p})) \equiv 1 - (1 - \tilde{p})^N$$
   This prevents severe catastrophic cancellation at $\tilde{p} \sim 10^{-4}$ and $N \sim 10^4$.
2. **Singleton / Support Protections:**
   - When $f_1 = N$, the Good–Turing coverage $C = 0$ singularity is avoided using the standard $f_1 := N - 1$ substitution.
   - When $f_1 = f_2 = 0$, the CWJ coverage bracket resolves $(0/0) \to 1$.
   - Explicit zeros in sparse matrices are purged using `Matrix::drop0()` so that `diff(@p)` accurately reflects observed feature support ($K$).
   - Non-negativity and integer checks are enforced inside the chunked processing blocks.
3. **Unit Consistency:** The Miller–Madow correction explicitly divides by $\ln(2)$ to convert nats to bits, avoiding the common literature trap.
4. **Hex Lattice Geometry:** Six hexagonal offsets on the Visium grid are validated against an independent $O(n^2)$ naive Moran's $I$ implementation, properly preserving the normalization factor $n/S_0$ when holes or isolated spots are present.

### B. Identified Flaws & Wrong Choices

#### 1. Statistical Inexactitude in `diagnose_correction_routes.R` (Line 194)
In [diagnose_correction_routes.R:L194](file:///home/alexej/Projects/Tirocinio/diagnose_correction_routes.R#L194):
```r
Spearman_rho_RouteM_vs_RouteI = cor(rank(route_m), route_i_rank)
```
- **The Issue:** `route_i_rank` is assigned from `st_rank$residual_y`, which contains the continuous numeric residuals from `lm.fit(rank(design), rank(plugin))`. These residuals are real numbers, not integer ranks.
- `cor()` without `method = "spearman"` computes the Pearson correlation.
- Consequently, this computes the Pearson correlation between the ranked Route M and the *unranked* continuous residuals of Route I (a semi-rank correlation).
- **The Correction:** It should either be:
  ```r
  Spearman_rho_RouteM_vs_RouteI = cor(route_m, route_i_rank, method = "spearman")
  # or equivalently:
  Spearman_rho_RouteM_vs_RouteI = cor(rank(route_m), rank(route_i_rank))
  ```

#### 2. Factual Contradiction in `docs/DECISIONS.md` on Criterion 4
In [DECISIONS.md:L133](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md#L133), the text summarizes the 16-sample cohort results:
> *"With criterion 4 passing everywhere and criterion 2 passing everywhere, three of the five A0 criteria separate nothing."*

- **The Reality:** In [cohort_estimator_criteria_summary.csv:L11](file:///home/alexej/Projects/Tirocinio/results/statistical_tests/cohort_estimator_criteria_summary.csv#L11) and [cohort_estimator_criteria.csv](file:///home/alexej/Projects/Tirocinio/results/statistical_tests/cohort_estimator_criteria.csv), `sample9` for `chao_shen` yields:
  $$\text{IQR} = 0.09911 \text{ bits} < 0.10 \text{ bits}$$
  Therefore, `C4_Dynamic_Range` is **FALSE** on `sample9` (15/16 pass, 93.8%).
- While Criterion 4 passed on all 3 initial diagnostic samples (`sample1`, `sample4`, `sample21`), stating that it passed *everywhere* across the 16-sample cohort is factually false.

#### 3. Internal Inconsistency in `docs/ESTIMATOR_PIVOT_PLAN.md`
- **Line 7:** *"Status: A0, A1, A2, B1, A3 and A4 are done; A5 (the estimator decision) is the next step and is deliberately left open... Everything from A3 onward is still forward-looking."*
- **Line 179:** *"### A5 — Gate and commit ✅ DONE"*
- The document header was left unedited after Stage A5 was completed and written up in the body.

#### 4. Pipeline Ordering Violation (Premature Stage C Execution)
- Section "Ordering" in [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md) specifies:
  $$\text{A5} + \text{B4} \longrightarrow \text{Stage C}$$
- Stage C explicitly states:
  > *"Do this once, after both A5 and B4 have landed. The whole point of doing A and B on three diagnostic samples is to avoid two 16-sample runs."*
- Despite this rule, [analyze_entropy.R](file:///home/alexej/Projects/Tirocinio/analyze_entropy.R) and [analyze_stemness.R](file:///home/alexej/Projects/Tirocinio/analyze_stemness.R) were modified to include `entropy_chao_shen` and executed on `sample1` (which fit SpaNorm, taking substantial memory and runtime).
- If Stage B concludes that SpaNorm should be dropped in favor of `NormalizeData()`, the edits to `analyze_entropy.R` will need to be re-done, and the `sample1` SpaNorm fit was unnecessary.

#### 5. QC Metric Schema Drift
- In [analyze_entropy.R:L191-197](file:///home/alexej/Projects/Tirocinio/analyze_entropy.R#L191-L197), seven new columns (`*_Chao_Shen_Entropy`) were added to `<sample>_qc_metrics.csv` while keeping the original seven `*_Raw_Plugin_Entropy` columns.
- However, [DECISIONS.md D1](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md#L15) still documents only the seven `*_Raw_Plugin_Entropy` columns.
- The plan indicated that columns should be renamed at Stage C; if both are to be retained permanently, D1 must be updated to document the new schema.

---

## 5. Verification & Diagnostic Cross-Check

Independent execution of the unit test suites confirms:
1. **Entropy Estimator Test Suite:**
   ```bash
   Rscript tests/test_entropy_estimators.R
   # Result: 49 passed, 0 failed
   ```
2. **Spatial Neighborhood Test Suite:**
   ```bash
   Rscript tests/test_spatial_neighbors.R
   # Result: 33 passed, 0 failed
   ```
3. **Cohort Statistical Outputs:**
   All 16 cohort samples have complete outputs for:
   - `<sample>_estimator_comparison.csv`
   - `<sample>_subsampling_stability.csv`
   - `<sample>_subsampling_reliability.csv`
   - `<sample>_correction_routes.csv`
   - Combined cohort summary CSVs in `results/statistical_tests/cohort_*`.

---

## 6. Recommendations for Next Actions

1. **Commit Working Tree Code (High Priority):**
   Stage and commit all untracked files ([R/load_sample.R](file:///home/alexej/Projects/Tirocinio/R/load_sample.R), [R/spatial_neighbors.R](file:///home/alexej/Projects/Tirocinio/R/spatial_neighbors.R), [R/diagnostic_sweep.R](file:///home/alexej/Projects/Tirocinio/R/diagnostic_sweep.R), `diagnose_*.R`, `tests/*`) to align git history with the commit logs.
2. **Synchronize Documentation Discrepancies:**
   - Update [ESTIMATOR_PIVOT_PLAN.md](file:///home/alexej/Projects/Tirocinio/docs/ESTIMATOR_PIVOT_PLAN.md) line 7 to reflect that A5 is completed.
   - Update [DECISIONS.md](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md) line 133 to state that C4 passed on 15 of 16 cohort samples (noting `sample9` at $0.0991$ vs $0.10$).
   - Update [DECISIONS.md D1](file:///home/alexej/Projects/Tirocinio/docs/DECISIONS.md#L15) to reflect the additional Chao–Shen columns in `_qc_metrics.csv`.
3. **Fix Semi-Rank Correlation in `diagnose_correction_routes.R`:**
   Change line 194 to `cor(route_m, route_i_rank, method = "spearman")`.
4. **Execute Stage B Before Further Cohort Work:**
   Implement B2 (Moran's $I$ on covariates) and B3 (SpaNorm vs. `NormalizeData()` comparison on stemness scores) on `sample1`, `sample4`, and `sample21` to make the Stage B4 branch decision before running any cohort-wide pipelines.
