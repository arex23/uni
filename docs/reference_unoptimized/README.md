# Reference copies — pre-optimization

These are the **unmodified** versions of two pipeline files as they stood at
commit `a1117c0`, before the memory work described in the D1 and D2 "Memory
note" paragraphs of `docs/DECISIONS.md`. They are kept here as the plain,
easier-to-explain statement of the method. Nothing in the pipeline reads this
directory — these files are reference material only.

| File here | Live counterpart |
|---|---|
| `analyze_entropy.R` | `../../analyze_entropy.R` |
| `shannon_entropy.R` | `../../R/shannon_entropy.R` |

There is no reference copy of `R/spanorm_lowmem.R`, because it did not exist
before the optimization. It has no methodological counterpart: deleting it *is*
the unoptimized state.

## What the live versions changed

Three things, all implementation-only. The numerical output is bit-identical —
verified with `identical()` against these files on the real sample1 object
(`data` layer, `shannon_entropy`, `shannon_entropy_raw`; maximum absolute
difference 0).

1. `analyze_entropy.R` sources `R/spanorm_lowmem.R` and calls
   `enable_spanorm_lowmem()`, which swaps SpaNorm's `logpac` adjustment for a
   kernel that does the same arithmetic in blocks of genes.
2. `analyze_entropy.R` drops the `SpatialExperiment` reference before the
   dense-to-sparse conversion of `logcounts`, so both representations are never
   in memory at once.
3. `R/shannon_entropy.R` evaluates its sparse branch over blocks of columns
   instead of over the whole matrix in one pass.

Peak RSS on sample1: crashed at ~15 GB → **7.59 GB**.

## Restoring the unoptimized pipeline

```sh
cp docs/reference_unoptimized/analyze_entropy.R  analyze_entropy.R
cp docs/reference_unoptimized/shannon_entropy.R  R/shannon_entropy.R
rm R/spanorm_lowmem.R
```

Be aware that this is the version that gets killed by the OOM killer on a
15.1 GB machine: SpaNorm's adjustment step alone holds about ten dense
15,556 × 11,709 double matrices at once, roughly 13.6 GB.

## For the write-up

The Methods section describes SpaNorm's published algorithm, and that algorithm
is the same in both versions — the block loop changes only the order in which
memory is allocated, not what is computed. So there is no need to run the
unoptimized code to be able to describe it: the optimized pipeline produces the
identical object, and the difference is already recorded as an implementation
note under D1 and D2 rather than as a methodological choice.
