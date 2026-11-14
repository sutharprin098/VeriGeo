# VeriGeo Methodology and Code Correction Log

## 1. Algorithmic and Mathematical Improvements
- **Iterative Collinearity Pruning in GWR (`protocol.py`)**:
  - *Previous*: Single-pass pruning calculated maximum pairwise correlation once and dropped one column. When 3+ predictors were collinear (e.g. `x1`, `x2`, `x3`), condition number remained $\kappa > 30$ ($66,186.9$).
  - *Corrected*: Implemented an iterative while-loop. At each step, the design matrix condition number $\kappa$ is computed on the current subset of standardized predictors. If $\kappa > 30$, the predictor with the highest maximum correlation against remaining columns is pruned. This repeats until $\kappa \le 30.0$ or a single predictor remains ($\kappa = 1.00$).
- **Compactness-Aware Sliver Filtering (`protocol.py`)**:
  - *Previous*: Evaluated pure strip width $2A/P < 1.0\text{ m}$, which pruned both narrow slivers and small compact footprints (e.g., $1.5\text{ m} \times 1.5\text{ m}$ utility vaults with $2A/P = 0.75\text{ m}$).
  - *Corrected*: Required both $W_{\text{eq}} = 2A/P < 1.0\text{ m}$ AND Polsby-Popper compactness $\text{PP} < 0.25$, preserving compact vaults while removing elongated digitizing slivers.
- **Metric Projection Routing**:
  - *Corrected*: Conformal UTM zone selection strictly enforced for distance and buffer operations to prevent high-latitude directional distortion, with Polar Stereographic fallback (`EPSG:3413`/`3031`) for $|\phi| \ge 80^\circ$. Equal-area projections (`EPSG:5070`, `EPSG:3035`, `EPSG:6933`) reserved for area, aggregation, and dissolve.
- **Canonical Provenance Digest**:
  - *Corrected*: Clarified that SHA-256 provides canonical hash integrity and execution trace provenance, replacing overclaiming terms like "signed proof" with "Spatial Verification Record (SPC)".

## 2. Experimental and Benchmark Results
- **Full 70-Task Execution Re-Run**:
  - Total runs logged: $535$ ($3 \times 70$ main methods $+ 10$ area-weighted joins $+ 15$ compactness-rule runs $+ 5 \times 60$ ablations).
  - Unverified GeoPandas (`naive` / `naive_pre`): $10/70$ passed ($14.3\%$, $0/60$ on fault-injected tasks, $10/10$ clean controls).
  - VeriGeo (`ns_vgp`): $67/70$ passed ($95.7\%$, $95\%$ Wilson CI $[88.1\%, 98.5\%]$), passing $15/15$ CRS, $13/15$ MAUP, $14/15$ TOPO, $15/15$ STAT, and $10/10$ CONTROL.
  - End-to-end execution latency: $25.44\text{ ms}$ mean runtime across all tasks.
