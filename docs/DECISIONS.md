# Decisions log

This file is the methodology record for the project: what was decided,
why, what alternatives were considered, and what happened when the
implementation was reviewed. Treat it as a running draft of the thesis
Methods section — write entries in full sentences, not just keywords.

Update `Status` as things progress: `draft` → `implemented` →
`under review` → `finalized`.

---

## D1 — Normalization & preprocessing
**Script:** `analyze_entropy.R`
**Decision:** [What normalization is applied, in what order, to what
input (raw counts / spots / genes filtered how)]
**Rationale:** [Why this approach, in this order]
**Alternatives considered / rejected:** [What else was possible, why not]
**Status:** draft
**Review notes:**

---


## D3 — Entropy–stemness correlation test
**Script:** `entropy_correlation.R` (called from `analyze_entropy.R`)
**Decision:** [Which statistical test, comparing entropy against what
stemness measure, at what unit of observation]
**Rationale:** [Why this test is appropriate — note: spatial data often
violates independence assumptions; state explicitly how/whether spatial
autocorrelation is accounted for]
**Alternatives considered / rejected:**
**Status:** draft
**Review notes:**

---

## D2 — Shannon entropy formulation
**Script:** `shannon_entropy.R` (called from analyze_entropy.R)
**Decision:** [Entropy computed over what unit — per spot? per gene? —
using what distribution, what log base, over which gene set (all genes?
top variable genes? a fixed panel?)]
**Rationale:** [Why this formulation captures "stemness-relevant"
heterogeneity]
**Alternatives considered / rejected:**
**Status:** draft
**Review notes:**

---

## D4 — Stemness gene reference / scoring
**Script:** `analyze_stemness.R`
**Decision:** [What gene set or scoring method defines "stemness" —
e.g. a published signature, a computed score like CytoTRACE-style, etc.
— and how it's applied per spot/region]
**Rationale:** [Why this reference is a valid ground truth for stemness
in meningioma specifically]
**Alternatives considered / rejected:**
**Status:** draft
**Review notes:**

---

## Template for new entries

```
## D<N> — <short title>
**Script:** <file>
**Decision:**
**Rationale:**
**Alternatives considered / rejected:**
**Status:** draft
**Review notes:**
```
