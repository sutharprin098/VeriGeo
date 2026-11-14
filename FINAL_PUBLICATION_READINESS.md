# Publication Readiness Review — Submission 6

## Overall assessment

**Status: conditional preprint candidate; not ready to describe as journal-ready or “minor revision.”**

The latest repository benchmark records 535 method executions over 70 synthetic tasks. The current `ns_vgp` results are 67/70 passed (95.7%; 95% Wilson CI 88.1–98.5), consisting of 57/60 fault-injected tasks and 10/10 clean controls. Category results are CRS 15/15, MAUP 13/15, topology 14/15, and spatial-statistics pre-checks 15/15. Mean `ns_vgp` latency over all 70 tasks is 25.44 ms. Three ground-truth tasks remain failures: MAUP-14, MAUP-15, and TOPO-15.

These figures describe the supplied deterministic benchmark; they are not estimates of performance on general GIS workloads. The STAT category is explicitly a pre-check/output test, not a model-level validation of GWR, LISA, or Getis-Ord.

## What is now supported by the evidence

1. The CRS benchmark compares the implementation against WGS84 ellipsoidal geodesics and local AEQD reference operations.
2. The MAUP benchmark uses 40,000 synthetic points and compares zone estimates with counts and attributes independently aggregated from the same points.
3. The repository reports the task-level CSV, summary JSON, and generated figures; the 535 rows are multiple methods/ablations over 70 tasks, not 535 independent tasks.
4. Iterative predictor pruning now resolves the three-predictor STAT-13 design-matrix diagnostic by dropping two predictors and leaving one (`kappa = 1.00`).
5. The manuscript discloses the three unresolved benchmark failures and narrows the statistical claims to the operations actually implemented.

## Critical limitations that must remain visible

- **No GWR/LISA/Getis-Ord estimator fit:** the shared statistics executor emits global/local Moran-style diagnostics for the STAT task labels. It does not fit GWR, LISA, or Getis-Ord models.
- **No permutation inference:** for `10 <= N < 30`, the pre-check records a conditional-permutation mode and a requested count of 999. The current executor does not run permutations or compute permutation-based p-values.
- **Synthetic scope:** all 70 tasks are generated from parameterized templates under a fixed seed. No empirical census/cadastral dataset or live LLM-agent workflow is evaluated.
- **Evaluator alignment:** CRS and MAUP have external reference truth. The topology short-width criterion is closely aligned with the repair threshold, and the STAT scorer tests finite Moran-style outputs, islands, and matrix conditioning rather than model estimation.
- **Partial target coverage:** the areal-interpolation routine normalizes allocations over intersecting target pieces and does not first reject incomplete source-to-target coverage. The caller must ensure the target zones cover the intended source support.
- **CRS units:** projected CRSs with non-metre units are not normalized before values such as `distance_m` are interpreted.
- **Ambiguous numeric columns:** an opaque non-negative integer field can still be inferred as extensive; MAUP-15 demonstrates this.
- **Narrow legitimate features:** the current thinness/elongation filter still removes the 0.80 m utility canal in TOPO-15.

## Required work before a journal submission

1. Add empirical census-block/tract and cadastral datasets with documented provenance and licensing.
2. Implement actual GWR, LISA, and Getis-Ord estimators (or remove those operator claims and task labels), and implement genuine permutation inference if it is to be claimed.
3. Add explicit target-coverage diagnostics and reject/flag partial coverage before full source totals are redistributed.
4. Add CRS axis-unit normalization and fail-closed semantic typing for ambiguous integer-valued columns.
5. Make sliver decisions feature-aware or use an analyst-configurable preserve/ignore policy for narrow infrastructure.
6. Execute the newly added `test_stat13_iterative_condition_precheck_only()` regression test for the three-predictor STAT-13 path. It is present on the revision branch but has not yet been run.
7. Rebuild the Overleaf/PDF bundle from the revised branch and inspect every rendered page before submission.

## Venue guidance

An EarthArXiv preprint may be reasonable after the revised manuscript and bundle are rebuilt and the limitations above are stated plainly. A peer-reviewed journal submission needs materially stronger external validation and genuine model-level statistical implementation. The available evidence does not support a “minor revision” readiness claim or a numerical acceptance/readiness score.
