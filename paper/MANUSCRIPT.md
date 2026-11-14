# NS-VGP: A Neuro-Symbolic Protocol for Verifying CRS Projections, MAUP Aggregation, Topology, and Spatial Econometrics in Autonomous GIS

**Author:** Prince Suthar  
**Affiliation:** AutoGIS Research Lab  
**Protocol Version:** `NS-VGP/1.0.0` (`VeriGeo 1.0`)  
**Repository:** https://github.com/sutharprin098/VeriGeo  
**Validation Data:** `benchmarks/per_task_validation_log.csv` (180 task runs) and `benchmarks/benchmark_results.json`

---

## Abstract

Autonomous GIS agents such as LLM-Geo, GIS Copilot, and GeoGPT use large language models to write and run spatial analysis scripts. Service frameworks such as Geospatial Agentic Services (GAS 1.0.0) wrap these agents behind standardized HTTP endpoints. In both designs, an execution is marked successful whenever the generated Python script finishes with exit code `0` and writes an output file. Because geometric libraries like GeoPandas and Shapely allow operations that violate cartographic and statistical assumptions without raising runtime exceptions, existing agents routinely produce flawed maps and tables. Four silent failures occur repeatedly in practice:
1. Planar Euclidean distance and area operations are executed directly on unprojected geographic coordinates (`EPSG:4326`), which introduces latitude-dependent secant distortion.
2. Spatial joins and polygon dissolves sum (`SUM`) intensive ratio variables such as median income, population density, $\text{PM}_{2.5}$, or NDVI while duplicating extensive population counts across intersecting zones.
3. Polygon overlays leave self-intersecting rings and thin boundary slivers in output layers.
4. Spatial autocorrelation and regression models (Moran's $I$, LISA, and GWR) run on disconnected contiguity graphs, small samples ($N < 30$), or collinear design matrices.

We present the **Neuro-Symbolic Verifiable GeoAI Protocol (NS-VGP / VeriGeo 1.0)**, a runtime verification and repair layer for autonomous GIS workflows. NS-VGP checks four formal conditions before and after each spatial operator: the Metric and Geodesic CRS Invariant ($\mathcal{A}_{\text{CRS}}$), the Semantic Scale Type and Pycnophylactic Mass-Conservation Invariant ($\mathcal{A}_{\text{MAUP}}$), the Topological 9-Intersection and Polsby-Popper Manifold Invariant ($\mathcal{A}_{\text{Topo}}$), and the Spatial Statistical Pre-Condition Invariant ($\mathcal{A}_{\text{Stat}}$). When a pre-condition or post-condition fails, a deterministic repair routine reprojects layers to local UTM or Equal-Area coordinates (`EPSG:5070` or `EPSG:6933`), replaces naive spatial joins with Tobler's iterative pycnophylactic interpolation for extensive counts and area-weighted averaging for intensive ratios, repairs invalid rings, removes slivers with Polsby-Popper compactness below $0.005$, and replaces disconnected contiguity matrices with row-standardized metric $k$-nearest neighbor weights. Each execution outputs a JSON **Spatial Proof Certificate (SPC)** signed with a SHA-256 hash. We evaluated NS-VGP on $60$ deterministic workflows ($420$ total runs including ablations) on Python `3.11.9`, GeoPandas `1.1.4`, and NumPy `1.26.4`. Unconstrained LLM-Geo and standard GAS 1.0.0 baselines finished $50.0\%$ and $83.3\%$ of tasks without crashing, but both scored $0.0\%$ on formal cartographic validity due to a $48.25\%\text{--}53.25\%$ mean area error, an $86.914\%$ population duplication error, and a $4.00\times$ inflation of intensive attributes. NS-VGP satisfied all four invariants across all $60$ tasks ($100.0\%$ validity), reduced mass-transfer error to $0.000\%$, and added a mean verification time of $67.32\text{ ms}$.

---

## 1. Introduction

Over the past three years, researchers have connected large language models (LLMs) directly to spatial processing libraries to automate GIS workflows (Li & Ning, 2023; Janowicz et al., 2020; Mai et al., 2024). Systems such as LLM-Geo (Li & Ning, 2023), GeoGPT (Zhang et al., 2024), GIS Copilot (Akinboyewa et al., 2025), and handbook-driven data retrieval agents (Ning et al., 2025) translate natural-language instructions into executable GeoPandas, Rasterio, PySAL (Rey et al., 2022), or PyQGIS code. In 2026, Li et al. (2026) introduced Geospatial Agentic Services (GAS 1.0.0), a web protocol (`GetCapabilities`, `DescribeAgent`, `ExecuteTask`) that lets clients invoke remote GIS agents using standardized JSON requests.

These systems evaluate correctness at the syntax and process level. If a generated script runs without throwing a Python exception and writes a non-empty vector or raster file, the orchestrator accepts the result. In spatial analysis, however, code that runs cleanly often violates basic rules of cartography, areal aggregation, topology, or spatial statistics.

We audited failure modes across open-source LLM GIS agents and service wrappers and grouped the most common errors into four categories:

1. **Unprojected Euclidean Metrics ($\mathcal{A}_{\text{CRS}}$):** When asked to buffer features by $1,000\text{ m}$ or compute polygon area, LLM agents frequently call `gdf.buffer()` or `gdf.area` directly on geographic coordinates (`EPSG:4326`), or align a second layer to the first layer's geographic CRS (`gdf2.to_crs(gdf1.crs)`). Because longitude degrees converge toward the poles by $\cos(\phi)$, planar area calculations on unprojected data distort areas by $\sec(\phi) - 1$, which reaches $32.10\%$ at $\phi = 40.8^\circ\text{N}$ and $128.94\%$ at $\phi = 64.1^\circ\text{N}$ (Snyder, 1987).
2. **MAUP Aggregation Errors ($\mathcal{A}_{\text{MAUP}}$):** When transferring census or environmental attributes across mismatched polygon boundaries (such as census tracts to watersheds), LLM scripts typically run `gpd.sjoin(zones, tracts, predicate="intersects")` followed by `.groupby(zone_id).sum()`. This single line commits two distinct errors. First, any source tract that crosses a target boundary is summed into every target zone it touches, inflating extensive count totals (`population`, `households`) by $116.6\%\text{--}148.3\%$ (Goodchild & Lam, 1980; Tobler, 1979). Second, summing intensive ratio variables (`median_income`, `pop_density`, `ndvi`, `pm25`) commits the ecological fallacy (Openshaw, 1984; Kuhn, 2012; Scheider & Ballatore, 2017), multiplying average income or vegetation indices by the number of intersecting polygons.
3. **Invalid Rings and Overlay Slivers ($\mathcal{A}_{\text{Topo}}$):** Real-world boundary files often contain self-intersecting rings ("bowtie" polygons). Furthermore, intersecting two slightly offset boundary layers produces narrow sliver polygons along shared edges (Egenhofer & Franzosa, 1991; Polsby & Popper, 1991). If left unpruned, these slivers distort feature counts and choropleth class breaks.
4. **Violated Spatial Statistical Assumptions ($\mathcal{A}_{\text{Stat}}$):** Agents running Moran's $I$, Local Indicators of Spatial Association (LISA; Anselin, 1995), Getis-Ord $G_i^*$ (Getis & Ord, 1992), or Geographically Weighted Regression (GWR; Brunsdon et al., 1996) default to Queen contiguity weights (`weights="queen"`) and normal $p$-value approximations. When the input layer contains disconnected island polygons (zero neighbors), missing values (`NaN`), small sample sizes ($N < 30$), or collinear predictors (condition number $\kappa(\mathbf{X}) > 30$), the model either crashes with a singular matrix error or reports invalid $z$-scores.

![Figure 1: NS-VGP Protocol Architecture](figures/fig1_protocol_architecture.png)

Prompt instructions alone do not reliably prevent these errors because LLMs do not inspect coordinate bounds, polygon intersection graphs, Polsby-Popper compactness ratios, or matrix condition numbers before emitting code. To solve this problem deterministically, we built the **Neuro-Symbolic Verifiable GeoAI Protocol (NS-VGP / VeriGeo 1.0)** (Figure 1). NS-VGP separates high-level task planning (handled by the LLM) from formal cartographic and statistical verification (handled by a symbolic engine in Python; d'Avila Garcez & Lamb, 2023).

---

## 2. Related Work

### 2.1 Autonomous GIS and Agent Interoperability
Li & Ning (2023) defined five goals for autonomous GIS: self-generating, self-organizing, self-verifying, self-executing, and self-growing. Follow-up systems implemented code generation for spatial analysis (Akinboyewa et al., 2025), automated data retrieval (Ning et al., 2025), and conversational tool calling (Zhang et al., 2024). In 2026, Li et al. (2026) introduced Geospatial Agentic Services (GAS 1.0.0) to standardize how GIS agents expose capabilities and exchange datasets over HTTP. GAS 1.0.0 focuses on service contracts and leaves internal verification to individual agent implementations. In practice, current reference agents only check whether a Python script throws an exception. NS-VGP implements the "self-verifying" requirement from Li & Ning (2023) as a reusable middleware layer that checks mathematical invariants on input and output layers.

### 2.2 Semantic Typing in Geoprocessing Workflows
Kuhn (2012) proposed core concepts of spatial information to formalize geographic meaning independent of software storage formats. Building on Kuhn's framework, Scheider & Ballatore (2017), Scheider et al. (2020), and Meerlo et al. (2021) developed Core Concept Data Types (CCD) to prevent invalid tool combinations during automated workflow synthesis. Their ontology distinguishes *extensive* attributes (which accumulate with spatial support area, such as total population) from *intensive* attributes (which are independent of area, such as population density or temperature). While prior work used these types in logic-based workflow planners (such as APE), NS-VGP applies intensive and extensive typing at runtime to intercept and repair spatial joins and polygon dissolves generated by LLMs.

### 2.3 MAUP and Pycnophylactic Areal Interpolation
Openshaw (1984) showed that aggregating point or polygon data into different zonal boundaries alters summary statistics and regression coefficients, a phenomenon known as the Modifiable Areal Unit Problem (MAUP). When transferring data between incongruent polygon sets, Goodchild & Lam (1980) and Goodchild et al. (1993) formulated area-weighted and dasymetric areal interpolation. Earlier, Tobler (1979) introduced smooth pycnophylactic interpolation, which computes a smooth density surface $\rho(x,y)$ over a grid while preserving the exact initial count within each source polygon ($\iint_{A_s} \rho(x,y)\,dx\,dy = C_s$, $\rho(x,y) \ge 0$). NS-VGP combines Tobler's iterative Laplacian surface with polygon intersection weights so that extensive counts are conserved to floating-point precision while intensive ratios are integrated by overlap area.

---

## 3. Formal Specification of the Four Spatial Axioms

We represent a single geoprocessing step as $W = (\text{op}, \mathcal{L}_1, \mathcal{L}_2, \Theta)$, where $\text{op}$ is the operator name, $\mathcal{L}_1 = (\mathcal{G}_1, \mathcal{X}_1, \text{CRS}_1)$ is the primary vector layer with geometries $\mathcal{G}_1 = \{g_1, \dots, g_M\}$, attribute table $\mathcal{X}_1$, and coordinate reference system $\text{CRS}_1$, $\mathcal{L}_2$ is an optional secondary layer, and $\Theta$ is the dictionary of operator parameters.

### 3.1 Axiom 1: Metric and Geodesic CRS Invariant ($\mathcal{A}_{\text{CRS}}$)

Let $\mathcal{O}_{\text{metric}}$ be the set of operators that compute Euclidean distances, areas, spatial weights, or polygon overlays (`buffer`, `area`, `length`, `distance`, `nearest`, `areal_interpolation`, `spatial_join`, `dissolve`, `intersection`, `morans_i`, `lisa`, `gwr`), and let $\mathcal{O}_{\text{area}} \subset \mathcal{O}_{\text{metric}}$ be the subset requiring equal-area preservation (`area`, `areal_interpolation`, `spatial_join`, `dissolve`). A workflow step $W$ satisfies $\mathcal{A}_{\text{CRS}}$ if and only if:

$$\text{op} \in \mathcal{O}_{\text{metric}} \implies \text{unit}(\text{CRS}_1) \in \{\text{m}, \text{ft}\} \land \neg\text{is\_geographic}(\text{CRS}_1), \qquad \mathcal{L}_2 \neq \varnothing \implies \text{CRS}_1 = \text{CRS}_2$$

**Distortion Calculation and Repair ($\mathcal{R}_{\text{CRS}}$):** When $\mathcal{L}_1$ is stored in `EPSG:4326` with bounding-box center $(\lambda_c, \phi_c)$ and spans $(\Delta\lambda, \Delta\phi)$, planar calculations in degree coordinates stretch East-West distances by the secant factor $\sec(\phi_c) = 1/\cos(\phi_c)$, producing a first-order areal distortion of $\delta_A(\phi_c) = (\sec(\phi_c) - 1) \times 100\%$. If $\mathcal{A}_{\text{CRS}}$ fails, $\mathcal{R}_{\text{CRS}}$ reprojects both $\mathcal{L}_1$ and $\mathcal{L}_2$ to $\text{CRS}^*$:
- `EPSG:5070` (NAD83 Conus Albers Equal-Area) if $(\text{op} \in \mathcal{O}_{\text{area}} \lor \Delta\lambda > 6.5^\circ)$ and the bounding box lies in North America;
- `EPSG:3035` (ETRS89 LAEA Europe) if $(\text{op} \in \mathcal{O}_{\text{area}} \lor \Delta\lambda > 6.5^\circ)$ and the bounding box lies in Europe;
- `EPSG:6933` (WGS 84 EASE-Grid 2.0 Global Equal-Area) for all other equal-area or wide-span operations;
- Local conformal UTM zone `EPSG:(32600 + Z)` (Northern Hemisphere) or `EPSG:(32700 + Z)` (Southern Hemisphere), where $Z = \lfloor (\lambda_c + 180)/6 \rfloor + 1 \in \{1,\dots,60\}$, for local distance and buffer operations.

### 3.2 Axiom 2: Semantic Scale Type and MAUP Mass Conservation ($\mathcal{A}_{\text{MAUP}}$)

The classifier $\tau(X_k) \in \{\tau_{\text{ext}}, \tau_{\text{int}}, \tau_{\text{nom}}\}$ (`verigeo/protocol.py`, lines 115–242) assigns each column in $\mathcal{X}_1$ using lexical rules followed by value-range inspection:
- **Nominal ($\tau_{\text{nom}}$):** Identifiers matching `id`, `fid`, `geoid`, `fips`, `code`, `name`, `type`, `zone`, or string/object columns.
- **Intensive ($\tau_{\text{int}}$):** Scale-invariant ratios, densities, rates, temperatures, and spectral indices matching `rate`, `density`, `pct`, `ratio`, `index`, `mean`, `median`, `per_capita`, `ndvi`, `evi`, `lst`, `temp`, `aqi`, `pm25`, `income`, `poverty`, or fractional values bounded in $[-1, 1]$ or $[0, 100]$. Intensive patterns are evaluated before extensive patterns so compound column names such as `pop_density` or `crime_rate` are classified as $\tau_{\text{int}}$.
- **Extensive ($\tau_{\text{ext}}$):** Additive totals matching `pop`, `population`, `total`, `count`, `households`, `cases`, `emissions`, `tons`, `crimes`, `jobs`, or non-negative integer counts.

For any zonal aggregation or polygon-to-polygon transfer from source polygons $\{A_s\}_{s=1}^M$ to target polygons $\{B_t\}_{t=1}^K$ with aggregation operator $\bigoplus_k$:

$$\forall X_k \text{ with } \tau(X_k) = \tau_{\text{int}}: \quad \bigoplus_k \neq \texttt{SUM}, \qquad \forall C_k \text{ with } \tau(C_k) = \tau_{\text{ext}}: \quad \left| \sum_{s \in \mathcal{S}_{\text{cov}}} C_{s,k} - \sum_{t=1}^{K} \hat{C}_{t,k} \right| \le \varepsilon_{\text{mass}}$$

where $\mathcal{S}_{\text{cov}} = \{s : |A_s \cap \bigcup_{t=1}^K B_t| > 0\}$ and $\varepsilon_{\text{mass}} = 10^{-5} \sum_{s \in \mathcal{S}_{\text{cov}}} |C_{s,k}|$.

**Pycnophylactic and Area-Weighted Repair ($\mathcal{R}_{\text{MAUP}}$):** Both layers are projected to an equal-area CRS. For each extensive column $C_k$, a density surface $\rho_k^{(r)}(x,y) \ge 0$ is computed on a $40 \times 40$ grid over $R = 35$ iterations of discrete 4-neighbor Laplacian smoothing ($\omega = 0.8$):

$$\tilde{\rho}_{i,j}^{(r)} = \max\!\left(0,\, (1-\omega)\rho_{i,j}^{(r-1)} + \frac{\omega}{4}\left(\rho_{i+1,j}^{(r-1)} + \rho_{i-1,j}^{(r-1)} + \rho_{i,j+1}^{(r-1)} + \rho_{i,j-1}^{(r-1)}\right)\right)$$

followed at each step by zone-wise mass rescaling $\rho_{i,j}^{(r)} = \tilde{\rho}_{i,j}^{(r)} \cdot C_{s,k} / \sum_{(u,v) \in A_s} \tilde{\rho}_{u,v}^{(r)}$. Let $c_{s,t}$ be the centroid of intersection piece $A_s \cap B_t$ with area $a_{s,t} = |A_s \cap B_t|$ and weight $w_{s,t} = a_{s,t} \cdot \max(\rho_k^{(R)}(c_{s,t}), 10^{-6})$. The repaired target values are:

$$\hat{C}_{t,k} = \sum_{s=1}^{M} C_{s,k} \frac{w_{s,t}}{\sum_{t'=1}^{K} w_{s,t'}} \quad (\tau_{\text{ext}}), \qquad \hat{Z}_{t,k} = \frac{\sum_{s=1}^{M} a_{s,t} Z_{s,k}}{\sum_{s=1}^{M} a_{s,t}} \quad (\tau_{\text{int}})$$

Summing $\hat{C}_{t,k}$ over $t = 1,\dots,K$ yields $\sum_{t=1}^K \hat{C}_{t,k} = \sum_{s \in \mathcal{S}_{\text{cov}}} C_{s,k}$, conserving covered population mass to floating-point precision while keeping $\hat{Z}_{t,k} \in [\min_s Z_{s,k}, \max_s Z_{s,k}]$.

### 3.3 Axiom 3: Topological 9-Intersection and Sliver Invariant ($\mathcal{A}_{\text{Topo}}$)

Every geometry $g \in \mathcal{G}$ must be non-empty, valid under the OGC Simple Features / Egenhofer 9-intersection rules (`g.is_valid == True`), and for polygon overlay outputs (`intersection`, `difference`, `sym_difference`), every polygon must satisfy either the Polsby-Popper compactness threshold or the relative area threshold:

$$\text{PP}(g) = \frac{4\pi \cdot \text{Area}(g)}{\text{Perimeter}(g)^2} \ge 0.005 \quad \lor \quad \frac{\text{Area}(g)}{\overline{\text{Area}}(\mathcal{G})} \ge 10^{-4}$$

### 3.4 Axiom 4: Spatial Statistical Pre-Conditions ($\mathcal{A}_{\text{Stat}}$)

For spatial statistical operators $\text{op} \in \{\text{morans\_i}, \text{lisa}, \text{getis\_ord}, \text{gwr}\}$ on target vector $\mathbf{y} \in \mathbb{R}^N$, spatial weights matrix $\mathbf{W} \in \mathbb{R}^{N \times N}$, and predictor matrix $\mathbf{X} \in \mathbb{R}^{N \times p}$:
1. $\mathbf{y}$ contains no `NaN` or $\pm\infty$ entries and $\text{Var}(\mathbf{y}) > 10^{-12}$.
2. Every spatial unit has at least one neighbor: $\forall i \in \{1,\dots,N\},\, \sum_{j=1}^N w_{ij} > 0$.
3. Sample size satisfies $N \ge 30$ for normal approximation, or conditional Monte Carlo permutation ($P \ge 999$) is used when $10 \le N < 30$.
4. For multivariate regression (`gwr`), the standardized design matrix $\tilde{\mathbf{X}}$ has condition number $\kappa(\tilde{\mathbf{X}}) < 30$.

---

## 4. Experimental Methodology and Empirical Results

### 4.1 Test Environment and Benchmark Protocol
All experiments were run on Windows using Python `3.11.9`, GeoPandas `1.1.4`, NumPy `1.26.4`, Shapely `2.0+`, PyProj `3.6+`, and Pytest `9.1.1` with deterministic seed `np.random.default_rng(42)`. The benchmark script `benchmarks/run_experiments.py` constructs $60$ tasks ($15$ per category) and runs each task under three modes (`llm_geo`, `gas_baseline`, `ns_vgp`) plus four single-axiom ablations (`no_A_CRS`, `no_A_MAUP`, `no_A_TOPO`, `no_A_STAT`), totaling $420$ executions. All $180$ main runs are logged row-by-row in `benchmarks/per_task_validation_log.csv`.

![Figure 2: Benchmark Comparison](figures/fig2_benchmark_comparison.png)

| Task Category ($N=15$ each) | Execution Mode | No-Crash Rate (%) | Verified Valid (%) | Mean Area Error (%) | Extensive Mass Error (%) | Intensive Inflation ($\times$) | Mean Latency ($\text{ms} \pm \text{SD}$) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. Metric & Geodesic CRS** | LLM-Geo (Unconstrained) | 100.0% | 0.0% | 53.25% | 0.000% | $1.00\times$ | $0.68 \pm 0.34$ |
| | GAS 1.0.0 (Standard Service) | 100.0% | 0.0% | 48.25% | 0.000% | $1.00\times$ | $0.21 \pm 0.09$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$24.25 \pm 17.99$** |
| **2. MAUP & Areal Interpolation** | LLM-Geo (Unconstrained) | 100.0% | 0.0% | 28.40% | 86.914% | $4.00\times$ | $21.48 \pm 11.24$ |
| | GAS 1.0.0 (Standard Service) | 100.0% | 0.0% | 28.40% | 86.914% | $4.00\times$ | $12.22 \pm 1.34$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$174.99 \pm 103.26$** |
| **3. Topological Manifold & Slivers** | LLM-Geo (Unconstrained) | 0.0% | 0.0% | — | — | — | $29.47 \pm 3.39$ |
| | GAS 1.0.0 (Standard Service) | 100.0% | 0.0% | 0.00% | 0.000% | $1.00\times$ | $0.29 \pm 0.04$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$50.99 \pm 7.31$** |
| **4. Spatial Statistical Pre-Conditions** | LLM-Geo (Unconstrained) | 0.0% | 0.0% | — | — | — | $0.16 \pm 0.02$ |
| | GAS 1.0.0 (Standard Service) | 33.3% | 0.0% | — | — | — | $0.19 \pm 0.06$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$19.04 \pm 1.83$** |
| **All Categories ($N=60$)** | LLM-Geo (Unconstrained) | 50.0% | 0.0% | 40.83% | 86.914% | $4.00\times$ | $12.95 \pm 14.12$ |
| | GAS 1.0.0 (Standard Service) | 83.3% | 0.0% | 38.33% | 86.914% | $4.00\times$ | $3.23 \pm 5.23$ |
| | **VeriGeo / NS-VGP (Ours)** | **100.0%** | **100.0%** | **0.00%** | **0.000%** | **$1.00\times$** | **$67.32 \pm 82.29$** |

---

### 4.2 Granular Task-by-Task Validation Log (`per_task_validation_log.csv`)

The table below shows exact empirical records extracted directly from `benchmarks/per_task_validation_log.csv`:

| Task ID | Operator & Setup | Unverified Output (`llm_geo` / `gas_baseline`) | NS-VGP Verified Output & Symbolic Repair Trace | Error Before $\to$ After | Certificate ID | Latency (ms) |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: |
| `CRS-01` | `clip_and_area` ($\phi=12.0^\circ\text{N}$) | CRS mismatch (`EPSG:4326` vs `3857`); $\sec(12^\circ)=1.0223$ | Reprojected layers from (`EPSG:4326`, `EPSG:3857`) to `EPSG:32611` (UTM) | $17.23\% \to 0.00\%$ | `SPC-D789FCF50DC0` | $45.93$ |
| `CRS-06` | `area` ($\phi=37.8^\circ\text{N}$) | Unprojected `EPSG:4326`; $\sec(37.8^\circ)=1.2656$ | Reprojected to `EPSG:5070` (CONUS Albers Equal-Area) | $26.56\% \to 0.00\%$ | `SPC-14882EF45229` | $81.47$ |
| `CRS-07` | `clip_and_area` ($\phi=40.8^\circ\text{N}$) | Unprojected `EPSG:4326`; $\sec(40.8^\circ)=1.3210$ | Reprojected layers from (`EPSG:4326`, `EPSG:3857`) to `EPSG:32618` | $47.10\% \to 0.00\%$ | `SPC-350B660AD388` | $26.82$ |
| `CRS-12` | `area` ($\phi=55.7^\circ\text{N}$) | Unprojected `EPSG:4326`; $\sec(55.7^\circ)=1.7745$ | Reprojected to `EPSG:6933` (Global Equal-Area) | $77.45\% \to 0.00\%$ | `SPC-919A0E014E7D` | $11.11$ |
| `CRS-15` | `buffer` ($1000\text{ m}$, $\phi=64.1^\circ\text{N}$) | Unprojected `EPSG:4326`; $\sec(64.1^\circ)=2.2894$ | Reprojected to local UTM Zone `EPSG:32628` | $128.94\% \to 0.00\%$ | `SPC-803EB03E6513` | $13.64$ |
| `MAUP-01` | `spatial_join` ($16 \to 9$ zones) | True pop $=79,902$, `sjoin` sum $=198,394$ ($+148.30\%$); True inc $=\$78,179$, `sjoin` sum $=\$296,699$ ($3.80\times$) | Reprojected to `EPSG:5070`; pycnophylactic pop $=79,902.000$; area-weighted inc $=\$78,179.0$ | $148.30\% \to 0.000\%$ | `SPC-534A31C1152B` | $213.92$ |
| `MAUP-02` | `spatial_join` ($16 \to 9$ zones) | True pop $=85,776$, `sjoin` sum $=200,236$ ($+133.44\%$); True inc $=\$81,838$, `sjoin` sum $=\$326,194$ ($3.99\times$) | Reprojected to `EPSG:5070`; pycnophylactic pop $=85,776.000$; area-weighted inc $=\$81,838.0$ | $133.44\% \to 0.000\%$ | `SPC-4DA9C249FA8B` | $203.01$ |
| `MAUP-03` | `dissolve` ($16 \to 4$ districts) | Dissolve summed intensive `median_income`: $\$317,121$ vs true $\$79,280$ ($4.00\times$) | Blocked `SUM` on $\tau_{\text{int}}$; computed area-weighted district income $\$79,280.0$ | $4.00\times \to 1.00\times$ | `SPC-2499A81F1722` | $45.48$ |
| `MAUP-08` | `spatial_join` ($16 \to 9$ zones) | True pop $=91,149$, `sjoin` sum $=213,969$ ($+134.75\%$); True inc $=\$90,494$, `sjoin` sum $=\$372,158$ ($4.11\times$) | Reprojected to `EPSG:5070`; pycnophylactic pop $=91,149.000$; area-weighted inc $=\$90,494.0$ | $134.75\% \to 0.000\%$ | `SPC-13DC8248BDF8` | $169.95$ |
| `TOPO-01` | `intersection` ($4$ parcels) | Input bowtie self-intersection + boundary sliver ($\text{PP}=0.00008 < 0.005$) | Repaired 1 bowtie via `ST_MakeValid`; reprojected to `EPSG:32618`; pruned 1 sliver ($3 \to 2$ features) | Sliver $\to 0$ | `SPC-8598E46CE735` | $49.27$ |
| `STAT-01` | `gwr` ($N=16$, $2$ islands) | `NaN` at index 2; 2 Queen islands; collinear $\kappa(\tilde{\mathbf{X}})=120,146.1$ | Imputed median $52.422$; built $k$-NN ($k=4$) in `EPSG:32618`; pruned `x1` ($\kappa \to 1.00$); $999$ perms | Invalid $\to$ Valid | `SPC-104F4F4E330C` | $19.35$ |
| `STAT-06` | `morans_i` ($N=24$, $2$ islands) | `NaN` at index 2; 2 disconnected Queen contiguity islands | Imputed median $60.467$; built row-standardized $k$-NN ($k=5$) in `EPSG:32618`; $999$ perms | Invalid $\to$ Valid | `SPC-5C352FD41F5A` | $18.45$ |

---

### 4.3 Cartographic Case Study (`MAUP-01`) and Ablation Analysis

![Figure 3: MAUP Case Study](figures/fig3_maup_pycnophylactic_case_study.png)

Figure 3 shows the spatial distribution on task `MAUP-01`:
- **Panel (a):** The $16$ source census tracts have a true total population of $79,902$ and a mean median household income of $\$78,179$.
- **Panel (b):** Running `gpd.sjoin(zones, tracts, predicate="intersects")` followed by `.groupby("zone_id").sum()` counts every boundary-straddling tract multiple times, producing a total population of $198,394$ ($+148.30\%$ error).
- **Panel (c):** Summing `median_income` across intersecting tracts inflates the mean zone income to $\$296,699$ ($3.80\times$ inflation).
- **Panel (d):** NS-VGP reprojects both layers to `EPSG:5070`, runs $35$ iterations of Tobler's pycnophylactic smoothing for `population` and `households`, and computes area-weighted averages for `median_income`, `pop_density`, `pm25`, and `ndvi`. The resulting target layer has an exact total population of $79,902.0$ ($0.000\%$ error) and a mean income of $\$78,179$.

![Figure 4: Ablation and Latency](figures/fig4_ablation_and_latency.png)

Figure 4(a) plots the single-axiom ablation results across all $60$ workflows. Disabling $\mathcal{A}_{\text{CRS}}$ reduces verified validity from $100.0\%$ to $25.0\%$ because $45$ of the $60$ workflows (Categories 1, 2, and 4) depend on metric coordinates. Disabling $\mathcal{A}_{\text{MAUP}}$, $\mathcal{A}_{\text{Topo}}$, or $\mathcal{A}_{\text{Stat}}$ each drops validity from $100.0\%$ to $75.0\%$. Figure 4(b) shows that the verification and repair runtime is $24.25 \pm 17.99\text{ ms}$ for $\mathcal{A}_{\text{CRS}}$, $174.99 \pm 103.26\text{ ms}$ for $\mathcal{A}_{\text{MAUP}}$, $50.99 \pm 7.31\text{ ms}$ for $\mathcal{A}_{\text{Topo}}$, and $19.04 \pm 1.83\text{ ms}$ for $\mathcal{A}_{\text{Stat}}$ ($67.32\text{ ms}$ overall mean).

---

## 5. References

1. **Li, Z., Khosravi Kazazi, A., Liu, R., Akinboyewa, T., Ning, H., et al. (2026).** Geospatial Agentic Services: A Framework for Interoperable Geospatial Intelligence. *Annals of GIS*. DOI: [10.1080/19475683.2026.2738374](https://doi.org/10.1080/19475683.2026.2738374)
2. **Li, Z., & Ning, H. (2023).** Autonomous GIS: the next-generation AI-powered GIS. *International Journal of Digital Earth*, 16(2), 4668–4686. DOI: [10.1080/17538947.2023.2278895](https://doi.org/10.1080/17538947.2023.2278895)
3. **Akinboyewa, T., Li, Z., Ning, H., & Lessani, M. N. (2025).** GIS Copilot: towards an autonomous GIS agent for spatial analysis. *International Journal of Digital Earth*, 18(1). DOI: [10.1080/17538947.2025.2497489](https://doi.org/10.1080/17538947.2025.2497489)
4. **Ning, H., Li, Z., Akinboyewa, T., & Lessani, M. N. (2025).** An autonomous GIS agent framework for geospatial data retrieval. *International Journal of Digital Earth*, 18(1). DOI: [10.1080/17538947.2025.2458688](https://doi.org/10.1080/17538947.2025.2458688)
5. **Zhang, Y., Wei, C., He, Z., & Yu, W. (2024).** GeoGPT: An assistant for understanding and processing geospatial tasks. *International Journal of Applied Earth Observation and Geoinformation*, 131, 103976. DOI: [10.1016/j.jag.2024.103976](https://doi.org/10.1016/j.jag.2024.103976)
6. **Janowicz, K., Gao, S., McKenzie, G., Hu, Y., & Bhaduri, B. (2020).** GeoAI: spatially explicit artificial intelligence techniques for geographic knowledge discovery and beyond. *International Journal of Geographical Information Science*, 34(4), 625–636. DOI: [10.1080/13658816.2019.1684500](https://doi.org/10.1080/13658816.2019.1684500)
7. **Mai, G., Huang, W., Sun, J., Song, S., Mishra, D., et al. (2024).** On the Opportunities and Challenges of Foundation Models for GeoAI (Vision Paper). *ACM Transactions on Spatial Algorithms and Systems*, 10(2), Article 11, 1–46. DOI: [10.1145/3653070](https://doi.org/10.1145/3653070)
8. **Tobler, W. R. (1979).** Smooth Pycnophylactic Interpolation for Geographical Regions. *Journal of the American Statistical Association*, 74(367), 519–530. DOI: [10.1080/01621459.1979.10481647](https://doi.org/10.1080/01621459.1979.10481647)
9. **Openshaw, S. (1984).** *The Modifiable Areal Unit Problem*. CATMOG No. 38, Geo Books.
10. **Goodchild, M. F., & Lam, N. S.-N. (1980).** Areal interpolation: a variant of the traditional spatial problem. *Geo-Processing*, 1(3), 297–312.
11. **Kuhn, W. (2012).** Core concepts of spatial information for transdisciplinary research. *International Journal of Geographical Information Science*, 26(12), 2267–2276. DOI: [10.1080/13658816.2012.722637](https://doi.org/10.1080/13658816.2012.722637)
12. **Scheider, S., & Ballatore, A. (2017).** Semantic typing of linked geoprocessing workflows. *International Journal of Digital Earth*, 11(1), 113–138. DOI: [10.1080/17538947.2017.1305457](https://doi.org/10.1080/17538947.2017.1305457)
13. **Scheider, S., Meerlo, R., Kasalica, V., & Lamprecht, A.-L. (2020).** Ontology of core concept data types for answering geo-analytical questions. *Journal of Spatial Information Science*, 20, 167–201. DOI: [10.5311/JOSIS.2020.20.620](https://doi.org/10.5311/JOSIS.2020.20.620)
14. **Egenhofer, M. J., & Franzosa, R. D. (1991).** Point-set topological spatial relations. *International Journal of Geographical Information Systems*, 5(2), 161–174. DOI: [10.1080/02693799108927841](https://doi.org/10.1080/02693799108927841)
15. **Polsby, D. D., & Popper, R. D. (1991).** The Third Criterion: Compactness as a Procedural Safeguard Against Partisan Gerrymandering. *Yale Law & Policy Review*, 9(2), 301–353.
16. **Anselin, L. (1995).** Local Indicators of Spatial Association—LISA. *Geographical Analysis*, 27(2), 93–115. DOI: [10.1111/j.1538-4632.1995.tb00338.x](https://doi.org/10.1111/j.1538-4632.1995.tb00338.x)
17. **Brunsdon, C., Fotheringham, A. S., & Charlton, M. E. (1996).** Geographically Weighted Regression: A Method for Exploring Spatial Nonstationarity. *Geographical Analysis*, 28(4), 281–298. DOI: [10.1111/j.1538-4632.1996.tb00936.x](https://doi.org/10.1111/j.1538-4632.1996.tb00936.x)
18. **d'Avila Garcez, A., & Lamb, L. C. (2023).** Neurosymbolic AI: The 3rd Wave. *Artificial Intelligence Review*, 56, 12387–12406. DOI: [10.1007/s10462-023-10448-w](https://doi.org/10.1007/s10462-023-10448-w)
