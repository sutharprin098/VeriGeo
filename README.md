# VeriGeo (`NS-VGP 1.0`): Neuro-Symbolic Verifiable GeoAI Protocol

<p align="left">
  <img src="https://img.shields.io/badge/Protocol-NS--VGP%2F1.0.0-10b981?style=flat-square" alt="Protocol NS-VGP 1.0.0">
  <img src="https://img.shields.io/badge/Python-3.10%2B-2563eb?style=flat-square" alt="Python 3.10+">
  <img src="https://img.shields.io/badge/GAS_1.0.0-Compatible-0284c7?style=flat-square" alt="GAS 1.0.0 Compatible">
  <img src="https://img.shields.io/badge/Verified_Validity-100%25-059669?style=flat-square" alt="100% Verified Validity">
  <img src="https://img.shields.io/badge/License-MIT-f59e0b?style=flat-square" alt="License MIT">
</p>

**VeriGeo (`NS-VGP 1.0`)** is a formal **Neuro-Symbolic Verifiable GeoAI Protocol** and self-healing execution engine for **Autonomous GIS** and **Geospatial Agentic Services (GAS 1.0.0)**.

While existing LLM-based GIS agents (e.g., *LLM-Geo*, *GIS Copilot*, *GeoGPT*) and standard service wrappers (*GAS 1.0.0*) achieve high syntactic code execution rates, they frequently commit silent yet catastrophic cartographic and statistical errors:
1. **Geodesic Metric Distortion ($\mathcal{A}_{\text{CRS}}$):** Running planar Euclidean `buffer`, `area`, `distance`, or `overlay` operators on unprojected geographic coordinates (`EPSG:4326`), incurring $30\text{--}100\%$ latitude-dependent secant area distortion ($\delta_A(\phi) = \sec(\phi) - 1$).
2. **MAUP & Ecological Fallacy ($\mathcal{A}_{\text{MAUP}}$):** Summing (`SUM`) scale-invariant *intensive* variables (`median_income`, `pop_density`, `poverty_rate`, `ndvi`, `pm25`) during spatial joins or polygon dissolves, while duplicating *extensive* counts (`population`, `households`) across intersecting polygons.
3. **Topological Manifold & Sliver Corruption ($\mathcal{A}_{\text{Topo}}$):** Failing on self-intersecting "bowtie" polygons or propagating microscopic overlay sliver artifacts.
4. **Spatial Econometric Pre-Condition Failure ($\mathcal{A}_{\text{Stat}}$):** Running Moran's $I$, LISA, or GWR on disconnected "island" contiguity matrices ($\sum_j w_{ij} = 0$), small samples ($N < 30$), or multicollinear predictors ($\kappa(X) > 30$).

**VeriGeo (`NS-VGP 1.0`)** interposes a **4-Axiom Symbolic GIScience Verifier & Self-Healing Critic** between neural planning and spatial execution, emitting a cryptographically signed **Spatial Proof Certificate (`SPC`)** with a SHA-256 provenance digest for every output.

---

## Architecture

![VeriGeo NS-VGP 1.0 Architecture](paper/figures/fig1_protocol_architecture.png)

---

## The Four Formal Spatial Axioms

| Axiom ID | Axiom Name | Formal Invariant | Deterministic Symbolic Auto-Repair ($\mathcal{R}$) |
| :--- | :--- | :--- | :--- |
| **$\mathcal{A}_{\text{CRS}}$** | **Metric & Geodesic CRS Invariant** | $\text{op} \in \mathcal{O}_{\text{metric}} \implies \text{unit}(\text{CRS}) \in \{\text{m}, \text{ft}\} \land \text{CRS}_1 \equiv \text{CRS}_2$ | Bounds $\sec(\phi)$ distortion; auto-reprojects to optimal Conformal UTM (`EPSG:326xx`/`327xx`) or Equal-Area (`EPSG:5070`/`3035`/`6933`). |
| **$\mathcal{A}_{\text{MAUP}}$** | **Semantic Scale Type & Pycnophylactic Mass Conservation** | $\forall X_k \in \tau_{\text{int}}: \bigoplus_k \neq \texttt{SUM} \;\land\; \left\|\sum_s C_s - \sum_t \hat{C}_t\right\| \le \varepsilon_{\text{mass}}$ | Classifies columns into $\tau_{\text{ext}}, \tau_{\text{int}}, \tau_{\text{nom}}$; runs **Tobler (1979) Smooth Pycnophylactic Interpolation** for extensive counts and **Area-Weighted Integral Averaging** for intensive variables. |
| **$\mathcal{A}_{\text{Topo}}$** | **Topological 9-Intersection & Manifold Invariant** | $\forall g \in \mathcal{G}: \text{is\_valid}(g) \land \text{PP}(g) = \frac{4\pi A}{P^2} \ge \tau_{\text{PP}}$ | Applies constructive manifold repair (`ST_MakeValid`) and prunes overlay micro-slivers via the **Polsby–Popper Compactness Index**. |
| **$\mathcal{A}_{\text{Stat}}$** | **Spatial Statistical Pre-Condition Invariant** | $\|\mathbf{y}\|_{\text{NaN}} = 0 \land \sigma^2(\mathbf{y}) > 0 \land \min_i \sum_j w_{ij} > 0 \land \kappa(\tilde{\mathbf{X}}) < 30$ | Imputes `NaN` via spatial median, heals disconnected Queen/Rook islands with row-standardized metric $k$-NN weights, switches $N < 30$ to $999$-permutation Monte Carlo inference, and prunes collinear predictors. |

---

## Empirical Benchmark ($N=60$ Workflows)

![Benchmark Comparison](paper/figures/fig2_benchmark_comparison.png)

| Task Category ($N=15$ each) | GeoAI Paradigm | Syntactic Exec. (%) | Verified Scientific Validity (%) | Mean Area Error (%) | Extensive Mass Error (%) | Intensive Inflation ($\times$) | Mean Latency ($\text{ms} \pm \text{SD}$) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. Metric & Geodesic CRS** | LLM-Geo (Unconstrained) | 100.0% | 0.0% | 49.33% | 0.00% | $1.00\times$ | $0.23 \pm 0.11$ |
| | GAS 1.0.0 (Standard Service) | 100.0% | 0.0% | 44.33% | 0.00% | $1.00\times$ | $0.54 \pm 0.23$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.00%** | **$1.00\times$** | **$40.57 \pm 24.36$** |
| **2. MAUP & Areal Interpolation** | LLM-Geo (Unconstrained) | 100.0% | 0.0% | 28.40% | 86.91% | $4.00\times$ | $24.96 \pm 14.76$ |
| | GAS 1.0.0 (Standard Service) | 100.0% | 0.0% | 28.40% | 86.91% | $4.00\times$ | $26.93 \pm 9.90$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$292.02 \pm 175.35$** |
| **3. Topological Manifold & Slivers** | LLM-Geo (Unconstrained) | 0.0% | 0.0% | — | — | — | $56.88 \pm 17.36$ |
| | GAS 1.0.0 (Standard Service) | 100.0% | 0.0% | 0.00% | 0.00% | $1.00\times$ | $0.33 \pm 0.07$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.00%** | **$1.00\times$** | **$73.00 \pm 17.83$** |
| **4. Spatial Statistical Pre-Conditions** | LLM-Geo (Unconstrained) | 0.0% | 0.0% | — | — | — | $0.52 \pm 0.29$ |
| | GAS 1.0.0 (Standard Service) | 33.3% | 0.0% | — | — | — | $0.16 \pm 0.03$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.00%** | **$1.00\times$** | **$19.37 \pm 3.42$** |
| **Overall Benchmark ($N=60$)** | LLM-Geo (Unconstrained) | 50.0% | 0.0% | 38.87% | 86.91% | $4.00\times$ | $20.65 \pm 25.85$ |
| | GAS 1.0.0 (Standard Service) | 83.3% | 0.0% | 36.37% | 86.91% | $4.00\times$ | $6.99 \pm 12.53$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$106.24 \pm 140.67$** |

---

## Cartographic Case Study: Preventing Ecological Fallacy & Mass Duplication

![MAUP Case Study](paper/figures/fig3_maup_pycnophylactic_case_study.png)

## 4-Axiom Ablation Study & Verification Latency

![Ablation and Latency](paper/figures/fig4_ablation_and_latency.png)

---

## Quickstart

### 1. Installation

```bash
git clone https://github.com/sutharprin098/VeriGeo.git
cd VeriGeo
pip install -e .[dev]
```

### 2. Verify & Auto-Heal a Spatial Workflow in Python

```python
import geopandas as gpd
from shapely.geometry import box
from verigeo import NeuroSymbolicGeoVerifier

# Source census tracts (Extensive: population | Intensive: median_income)
tracts = gpd.GeoDataFrame(
    {
        "tract_id": ["T1", "T2", "T3", "T4"],
        "population": [2500.0, 4500.0, 3000.0, 5000.0],
        "median_income": [50000.0, 70000.0, 60000.0, 80000.0],
    },
    geometry=[
        box(-77.90, 40.75, -77.85, 40.80),
        box(-77.85, 40.75, -77.80, 40.80),
        box(-77.90, 40.80, -77.85, 40.85),
        box(-77.85, 40.80, -77.80, 40.85),
    ],
    crs="EPSG:4326",
)

# Incongruent target planning zones
zones = gpd.GeoDataFrame(
    {"zone_id": ["Z1", "Z2", "Z3"]},
    geometry=[
        box(-77.90, 40.75, -77.83, 40.82),
        box(-77.83, 40.75, -77.80, 40.82),
        box(-77.90, 40.82, -77.80, 40.85),
    ],
    crs="EPSG:4326",
)

verifier = NeuroSymbolicGeoVerifier(auto_heal=True)
verified_gdf, spc = verifier.verify_and_execute(
    operator="spatial_join",
    primary_gdf=tracts,
    secondary_gdf=zones,
    parameters={"aggfunc": "sum"},  # Illegal SUM on intensive median_income is automatically intercepted & healed!
)

print("Verifiable Guarantee:", spc.verifiable_guarantee)
print("Conserved Population Total:", verified_gdf["population"].sum())  # Exact: 15000.0
print("Repairs Applied:", spc.repairs_applied)
print("Provenance SHA-256:", spc.provenance_sha256)
```

### 3. Reproduce Benchmark & Publication Figures

```bash
python benchmarks/run_experiments.py
pytest -v
```

---

## Repository Structure

```text
VeriGeo/
├── verigeo/
│   ├── __init__.py          # Public API exports
│   ├── protocol.py          # Core 4-Axiom Verifier, Pycnophylactic Solver & SPC Engine
│   ├── gas_middleware.py    # GAS 1.0.0 Middleware & Merkle-Chained DAG Executor
│   └── cli.py               # Command-line interface (`verigeo verify`)
├── benchmarks/
│   ├── run_experiments.py   # 60-workflow quantitative benchmark & figure generator
│   └── benchmark_results.json
├── paper/
│   ├── main.tex             # Two-column LaTeX manuscript (IJGIS / Annals of GIS / arXiv)
│   ├── references.bib       # Verified BibTeX bibliography with DOIs
│   ├── MANUSCRIPT.md        # Full Markdown manuscript with equations & tables
│   └── figures/             # 300-DPI publication figures (PNG + PDF)
├── tests/
│   └── test_verigeo.py      # Comprehensive pytest test suite
├── pyproject.toml
└── LICENSE
```

---

## Citation

```bibtex
@article{suthar2026verigeo,
  author  = {Suthar, Prince},
  title   = {{NS-VGP}: A Neuro-Symbolic Verifiable {GeoAI} Protocol for Formally Guaranteed, Mass-Conserving, and Scale-Aware Autonomous {GIS}},
  journal = {arXiv preprint / Working Paper},
  year    = {2026},
  url     = {https://github.com/sutharprin098/VeriGeo}
}
```
