**EarthArXiv Preprint · Non-peer-reviewed preprint**

# VeriGeo: Runtime Verification for Projections, Areal Interpolation, Topology, and Spatial Statistics in Automated GIS

**Prince Suthar**  
*AutoGIS Research Lab*  
`https://github.com/sutharprin098/VeriGeo`  
*November 2026*

> **Manuscript source:** This file is a concise Markdown companion to `main.tex`. The LaTeX manuscript is the canonical submission file. Update this summary whenever benchmark metrics or method claims change.

---

## Abstract

Automated GIS scripts can complete without runtime exceptions while still producing spatially invalid outputs. We present **VeriGeo (NS-VGP 1.0)**, a Python runtime checking and repair layer for coordinate reference systems, areal interpolation, polygon topology, and selected spatial-statistics preconditions. We evaluate it on a deterministic synthetic benchmark of 70 tasks (60 fault-injected tasks and 10 clean controls), with 535 logged method/ablation executions. The unverified GeoPandas baseline passed 0/60 fault-injected tasks and 10/10 controls (14.3% overall). VeriGeo passed 57/60 fault-injected tasks and 10/10 controls (67/70 = 95.7% overall; 95% Wilson CI 88.1–98.5%), with 25.44 ms mean runtime across all 70 tasks. Three tasks remain failures: MAUP-14 (unmasked water barrier), MAUP-15 (opaque integer-valued income attribute), and TOPO-15 (legitimate 0.80 m canal removed by the sliver filter). The STAT benchmark evaluates finite Moran-style outputs and selected preconditions; it does not validate fitted GWR, LISA, or Getis-Ord estimators or permutation-based inference.

## 1. Problem and contribution

A GIS script can exit with code zero while using an unsuitable CRS, double-counting extensive attributes across polygon boundaries, summing intensive attributes, retaining invalid/elongated overlay fragments, or running spatial-statistics routines with missing values and disconnected weights. VeriGeo wraps selected spatial operators in pre-check, repair, execution, and post-check steps and returns a structured **Spatial Verification Record (SPC)**.

The current SPC SHA-256 digest covers audit metadata (operator, stringified parameters, axiom results, repairs, and guarantee flag). It is not a hash of the complete input/output spatial datasets, a digital signature, or proof of scientific correctness.

## 2. Implemented checks and boundaries

| Check | Current behaviour | Important boundary |
| :--- | :--- | :--- |
| **CRS** | Selects a projected CRS for metric operations; uses equal-area projections for area/aggregation and UTM/polar options for distance/buffer operations. | Projected CRSs with non-metre axis units are not normalized before metre-valued parameters such as `distance_m` are interpreted. |
| **Attribute scale / MAUP** | Uses lexical and value-range heuristics to classify columns as nominal, intensive, or extensive. Uses pycnophylactic allocation for extensive values and area-weighted averages for intensive values. | Ambiguous non-negative integer fields can be classified as extensive. Area-weighted averaging is not generally an exact method for medians or socioeconomic rates. Full target coverage is not checked by the interpolation routine; callers must validate support coverage before use. |
| **Geometry / slivers** | Repairs invalid geometries with `make_valid()`. Removes overlay pieces only when equivalent width is below 1.0 m and Polsby-Popper compactness is below 0.25. | TOPO-15 shows that the current heuristic still removes a legitimate 0.80 m canal. It is not a complete implementation of the 9-intersection spatial-relation model. |
| **Spatial-statistics pre-checks** | Checks target finite values and variance, detects islands, applies a k-NN weights fallback, and iteratively prunes collinear predictors. | For 10 ≤ N < 30, the code records a conditional-permutation mode and requests 999 permutations, but does not run them or calculate p-values. The shared executor emits Moran-style diagnostics for GWR/LISA-labelled tasks; it does not fit GWR, LISA, or Getis-Ord models. |

## 3. Benchmark results

The benchmark uses WGS84 ellipsoidal geodesics and local AEQD reference operations for CRS tasks, 40,000 synthetic resident points for MAUP tasks, minimum rotated bounding rectangle widths for topology tasks, and finite Moran-style outputs / island / condition-number checks for STAT tasks. All 70 tasks are generated from parameterized templates under a fixed seed. The topology and STAT evaluators are closely aligned with parts of the implementation; these results should not be read as independent model-level validation.

| Category | Baseline result | VeriGeo result | VeriGeo mean latency |
| :--- | :--- | :--- | ---: |
| CRS distortion (15) | 0/15 passed; 35.00% mean geodesic area error | 15/15 passed; 0.024% mean geodesic area error | 11.61 ± 8.54 ms |
| MAUP aggregation (15) | 0/15 passed; 169.56% mean population overcount on joins | 13/15 passed; zero total population error; 10.10% mean zone NMAE across 10 joins (3.84% across nine Gaussian-hotspot joins) | 66.82 ± 32.54 ms |
| Topology / slivers (15) | 0/15 passed; 1.87 slivers per task on average | 14/15 passed; TOPO-15 remains a failure | 29.09 ± 3.90 ms |
| Spatial-statistics pre-checks (15) | 0/15 passed | 15/15 passed under finite-output, zero-island, and design-matrix checks only | 8.41 ± 1.80 ms |
| Clean controls (10) | 10/10 passed | 10/10 passed with zero false-positive repairs | 4.18 ± 0.78 ms |
| **All tasks (70)** | **10/70 passed (14.3%)** | **67/70 passed (95.7%)**; 57/60 fault-injected and 10/10 clean controls | **25.44 ms mean** |

## 4. Remaining failures and repaired stress case

1. **MAUP-14 — unmasked coastal water barrier:** Pycnophylactic smoothing diffuses population across an uninhabited water boundary. Total population remains conserved, but zone NMAE is 66.46%, compared with 56.27% for area-weighted interpolation.
2. **MAUP-15 — opaque `attr_val` column:** Integer-valued median income is inferred as an extensive count and summed, producing approximately 3.893× inflation. Explicit semantic overrides are needed for ambiguous columns.
3. **TOPO-15 — narrow infrastructure:** The 0.80 m utility canal is removed by the thinness/elongation filter; only 67.57% of legitimate area is retained.
4. **STAT-13 — three collinear predictors:** The revised iterative pre-check removes two predictors and leaves one, recording a condition number of 1.00. This is a pre-check result, not a successful GWR fit. A dedicated regression test has been added on the revision branch but has not yet been executed.

## 5. Reproducibility and next steps

The repository contains the benchmark generator, per-task CSV, JSON summary, test suite, and figure generation scripts. The previously recorded Windows test run reports 8 passed; STAT-13 and missing-CRS regression tests have since been added but both are pending execution. No fresh pytest run or LaTeX compilation was performed as part of this review.

Before journal submission, the next priorities are empirical census/cadastral validation, actual GWR/LISA/Getis-Ord implementations or narrowed operator claims, genuine permutation inference if claimed, partial-coverage validation, CRS axis-unit normalization, fail-closed semantic typing for ambiguous attributes, and infrastructure-aware sliver handling.
