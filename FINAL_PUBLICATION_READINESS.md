# Publication Readiness and Peer Review Assessment

## 1. Overall Readiness Evaluation

| Venue | Readiness Score (out of 10) | Recommended Verdict | Key Strengths & Remaining Action Items |
| :--- | :---: | :---: | :--- |
| **EarthArXiv Preprint** | **9.5 / 10** | **Ready for Submission** | Fully reproducible benchmark ($535$ logged runs), honest failure case documentation (`MAUP-14`, `MAUP-15`, `TOPO-15`), mathematically sound, verified citations with DOIs. |
| **IJGIS / Computers & Geosciences** | **8.5 / 10** | **Minor Revision** | High technical rigor; remaining item for full journal submission is expanding validation to empirical census block-to-tract data alongside synthetic benchmarks. |

## 2. Scientific Strengths
1. **Zero Circularity in Primary Benchmarks**: WGS84 geodesic area calculations via `pyproj.Geod` and Azimuthal Equidistant feature-centered projections provide ground truth completely external to VeriGeo.
2. **Explicit Disclosure of Real Failure Cases**:
   - `MAUP-14`: Shows where unmasked pycnophylactic smoothing diffuses population across water barriers, underscoring the necessity of dasymetric masks.
   - `MAUP-15`: Demonstrates the boundary condition of lexical heuristics on opaque column names (`attr_val`) without unit metadata.
   - `TOPO-15`: Demonstrates the trade-off in fixed geometric sliver thresholds on narrow infrastructure features ($0.80\text{ m}$ utility canals).
3. **Rigorous Statistical Accounting**: Exactly $535$ runs accounted for ($3 \times 70 + 10 + 15 + 5 \times 60$). Wilson score 95% confidence intervals recomputed and exact to stated decimals.
4. **Iterative Multicollinearity Resolution**: Fixed single-pass collinearity pruning by implementing iterative eigenvalue condition number diagnostics ($\kappa \le 30.0$).
5. **Verified Foundational Bibliography**: 13 foundational citations integrated with verified DOIs across spatial analysis, map projections, change-of-support, and spatial econometrics.

## 3. Disclosed Limitations
- **Synthetic Task Generative Models**: The 70 benchmark tasks use synthetic geometric grids and point clusters to allow 100% deterministic reproducibility without paid API keys or commercial GIS licenses.
- **Fail-Open Heuristic Limits**: Attribute classification based on lexical tokens requires user confirmation or explicit schema definitions when working with legacy datasets with uninformative attribute headers.
- **Unmasked Discontinuity Smoothing**: Pycnophylactic smoothing assumes spatial continuity across tract borders and should be paired with land-use masks when natural barriers divide administrative units.
