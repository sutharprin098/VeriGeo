# NS-VGP: A Neuro-Symbolic Protocol for Verifying CRS Projections, MAUP Aggregation, Topology, and Spatial Econometrics in Autonomous GIS

**Author:** Prince Suthar  
**Affiliation:** AutoGIS Research Lab  
**Protocol Version:** `NS-VGP/1.0.0` (`VeriGeo 1.0`)  
**Repository:** https://github.com/sutharprin098/VeriGeo  
**Validation Data:** `benchmarks/per_task_validation_log.csv` (535 task-method executions across 70 tasks) and `benchmarks/benchmark_results.json`

---

## Abstract

Autonomous GIS agents such as LLM-Geo, GIS Copilot, and GeoGPT use large language models to write and run spatial analysis scripts. In these autonomous agent architectures, an execution is marked successful whenever the generated Python script finishes with exit code `0` and writes an output file. Because geometric libraries like GeoPandas and Shapely allow operations that violate cartographic and statistical assumptions without raising runtime exceptions, generated scripts frequently produce flawed maps and tables. Four silent failures occur repeatedly in automated workflows:
1. Planar Euclidean distance and area operations are executed on unprojected geographic coordinates (`EPSG:4326`) or converted via fixed degree-to-meter constants (`111,320 m/deg`), introducing latitude-dependent secant distortion.
2. Spatial joins and polygon dissolves sum (`SUM`) intensive ratio variables such as median income, population density, $\text{PM}_{2.5}$, or NDVI while duplicating extensive population counts across intersecting zones.
3. Polygon overlays leave self-intersecting rings and narrow boundary slivers in output layers.
4. Spatial autocorrelation and regression models (Moran's $I$, LISA, and GWR) run on disconnected contiguity graphs, missing values, small samples ($N < 30$), or collinear design matrices.

We present the **Neuro-Symbolic Verifiable GeoAI Protocol (NS-VGP / VeriGeo 1.0)**, a runtime verification and repair layer for autonomous GIS workflows. NS-VGP checks four formal conditions before and after each spatial operator: the Metric and Geodesic CRS Invariant ($\mathcal{A}_{\text{CRS}}$), the Semantic Scale Type and Pycnophylactic Mass-Conservation Invariant ($\mathcal{A}_{\text{MAUP}}$), the Topological 9-Intersection and Mean-Width Manifold Invariant ($\mathcal{A}_{\text{Topo}}$), and the Spatial Statistical Pre-Condition Invariant ($\mathcal{A}_{\text{Stat}}$). When a pre-condition or post-condition fails, a deterministic repair routine reprojects layers to local UTM or Equal-Area coordinates (`EPSG:5070`, `EPSG:3035`, or `EPSG:6933`), replaces naive spatial joins with Tobler's iterative pycnophylactic interpolation for extensive counts and area-weighted averaging for intensive ratios, repairs invalid rings, removes boundary slivers whose mean width $2A/P$ falls below $1.0\text{ m}$, and replaces disconnected contiguity matrices with row-standardized metric $k$-nearest neighbor weights. Each execution outputs a JSON **Spatial Proof Certificate (SPC)** signed with a SHA-256 hash.

We evaluated NS-VGP against idiomatic GeoPandas baselines (`naive`, `naive_pre`, and standard area-weighted interpolation `areaw`) and single-axiom ablations across $70$ tasks ($60$ adversarial tasks plus $10$ clean projected control tasks, totaling $535$ logged runs) using an independent ground-truth evaluator based on WGS84 ellipsoidal geodesics (`pyproj.Geod`), $40,000$ synthetic individual point populations, and minimum rotated bounding rectangle widths. Unverified `naive` and `naive_pre` baselines passed $0/60$ ($0.0\%$) of adversarial tasks and $10/10$ ($100.0\%$) of clean control tasks ($14.3\%$ overall), exhibiting a $35.00\%$ mean ($26.87\%$ median) geodesic area error, a $166.05\%$ population duplication error on spatial joins, and a $4.00\times$ inflation of intensive attributes. On spatial joins, NS-VGP's pycnophylactic solver reduced per-zone population normalized mean absolute error (NMAE) against the underlying $40,000$-point ground truth to $3.62\%$, outperforming standard area-weighted interpolation ($7.92\%$ NMAE) by more than $2\times$ while conserving total population to floating-point precision ($1.93 \times 10^{-15}\%$ error). Across all $70$ tasks, NS-VGP passed $70/70$ ($100.0\%$, $95\%$ Wilson CI $[94.8\%, 100.0\%]$) with zero false-positive repairs on the $10$ control tasks and a mean latency of $68.16\text{ ms}$.

---

## 1. Introduction

Over the past three years, researchers have connected large language models (LLMs) to spatial processing libraries to automate GIS workflows (Li & Ning, 2023; Janowicz et al., 2020; Mai et al., 2024). Systems such as LLM-Geo (Li & Ning, 2023), GeoGPT (Zhang et al., 2024), GIS Copilot (Akinboyewa et al., 2025), and handbook-driven data retrieval agents (Ning et al., 2025) translate natural-language instructions into executable GeoPandas, Rasterio, PySAL (Rey et al., 2022), or PyQGIS code.

These systems evaluate correctness at the syntax and process level. If a generated script runs without throwing a Python exception and writes a non-empty vector or raster file, the orchestrator accepts the result. In spatial analysis, however, code that runs cleanly often violates basic rules of cartography, areal aggregation, topology, or spatial statistics.

We audited common GeoPandas patterns emitted by LLM coding assistants and grouped the most frequent silent errors into four categories:

1. **Unprojected Euclidean Metrics ($\mathcal{A}_{\text{CRS}}$):** When asked to buffer features by $1,000\text{ m}$ or compute polygon area on geographic coordinates (`EPSG:4326`), generated scripts often convert meters to degrees using the equatorial constant $1^\circ \approx 111,320\text{ m}$ (`gdf.buffer(1000 / 111320.0)`) or multiply square degrees by $111,320^2$. Because longitude degrees converge toward the poles by $\cos(\phi)$, this planar degree approximation distorts areas relative to the WGS84 ellipsoid by $6.05\%$ at $\phi = 18.2^\circ\text{N}$, $26.87\%$ at $\phi = 37.8^\circ\text{N}$, and $112.63\%$ at $\phi = 60.2^\circ\text{N}$ (Snyder, 1987).
2. **MAUP Aggregation Errors ($\mathcal{A}_{\text{MAUP}}$):** When transferring census or environmental attributes across incongruent polygon boundaries (such as census tracts to planning zones), standard scripts run `gpd.sjoin(zones, tracts, predicate="intersects")` followed by `.groupby(zone_id).sum()`. This single pattern commits two distinct errors. First, any source tract that crosses a target boundary is summed into every target zone it touches, inflating extensive count totals (`population`, `households`) by $152.95\%\text{--}183.04\%$ across our $10$ spatial join tasks (Goodchild & Lam, 1980; Tobler, 1979). Second, summing intensive ratio variables (`median_income`, `pop_density`, `ndvi`, `pm25`) during joins or dissolves commits the ecological fallacy (Openshaw, 1984; Kuhn, 2012; Scheider & Ballatore, 2017), multiplying average income by $3.89\times\text{--}4.32\times$.
3. **Invalid Rings and Overlay Slivers ($\mathcal{A}_{\text{Topo}}$):** Real-world boundary files frequently contain self-intersecting rings ("bowtie" polygons). Furthermore, intersecting two slightly offset boundary layers produces narrow sliver polygons along shared edges (Egenhofer & Franzosa, 1991; Polsby & Popper, 1991). While `.make_valid()` fixes self-intersections, it does not remove geometric slivers, which distort feature counts and choropleth class breaks.
4. **Violated Spatial Statistical Assumptions ($\mathcal{A}_{\text{Stat}}$):** Scripts running Moran's $I$, Local Indicators of Spatial Association (LISA; Anselin, 1995), Getis-Ord $G_i^*$ (Getis & Ord, 1992), or Geographically Weighted Regression (GWR; Brunsdon et al., 1996) typically default to Queen contiguity weights (`weights="queen"`). When the input layer contains disconnected island polygons (zero touching neighbors), missing values (`NaN`), small sample sizes ($N < 30$), or collinear predictors (condition number $\kappa(\mathbf{X}) > 30$), naive scripts either propagate `NaN` values, leave zero-weight rows in the contiguity matrix, or invert ill-conditioned design matrices ($\kappa \in [61706, 133728]$).

![Figure 1: NS-VGP Protocol Architecture](figures/fig1_protocol_architecture.png)

Prompt instructions alone do not reliably prevent these errors because LLMs do not inspect coordinate bounds, polygon intersection graphs, sliver widths, or matrix condition numbers before emitting code. To address this deterministically, we built the **Neuro-Symbolic Verifiable GeoAI Protocol (NS-VGP / VeriGeo 1.0)** (Figure 1). NS-VGP separates high-level task planning from formal cartographic and statistical verification (d'Avila Garcez & Lamb, 2023).

---

## 2. Related Work

### 2.1 Autonomous GIS and Self-Verifying Agents
Li & Ning (2023) defined five goals for autonomous GIS: self-generating, self-organizing, self-verifying, self-executing, and self-growing. Follow-up systems implemented code generation for spatial analysis (Akinboyewa et al., 2025), automated data retrieval (Ning et al., 2025), and conversational tool calling (Zhang et al., 2024). Existing autonomous GIS implementations focus on code generation and leave runtime verification to Python exception checks. NS-VGP implements the "self-verifying" requirement from Li & Ning (2023) as a reusable Python verifier and DAG middleware layer that checks mathematical invariants on input and output layers.

### 2.2 Semantic Typing in Geoprocessing Workflows
Kuhn (2012) proposed core concepts of spatial information to formalize geographic meaning independent of software storage formats. Building on Kuhn's framework, Scheider & Ballatore (2017), Scheider et al. (2020), and Meerlo et al. (2021) developed Core Concept Data Types (CCD) to prevent invalid tool combinations during automated workflow synthesis. Their ontology distinguishes *extensive* attributes (which accumulate with spatial support area, such as total population) from *intensive* attributes (which are independent of area, such as population density or median income). While prior work used these types in logic-based workflow planners (such as APE), NS-VGP applies intensive and extensive typing at runtime to intercept and repair spatial joins and polygon dissolves.

### 2.3 MAUP and Pycnophylactic Areal Interpolation
Openshaw (1984) showed that aggregating point or polygon data into different zonal boundaries alters summary statistics and regression coefficients, a phenomenon known as the Modifiable Areal Unit Problem (MAUP). When transferring data between incongruent polygon sets, Goodchild & Lam (1980) and Goodchild et al. (1993) formulated area-weighted areal interpolation, which assumes uniform population density within each source tract. Earlier, Tobler (1979) introduced smooth pycnophylactic interpolation, which computes a smooth density surface $\rho(x,y)$ over a grid while preserving the exact initial count within each source polygon ($\iint_{A_s} \rho(x,y)\,dx\,dy = C_s$, $\rho(x,y) \ge 0$). When real populations cluster in spatial hotspots rather than uniform blocks, Tobler's smooth density surface captures cross-boundary gradients more accurately than uniform area weighting. NS-VGP combines Tobler's iterative Laplacian surface with exact polygon intersection weights so that extensive counts are conserved to floating-point precision while reducing zone-level estimation error relative to uniform area weighting.

---

## 3. Formal Specification of the Four Spatial Axioms

We represent a single geoprocessing step as $W = (\text{op}, \mathcal{L}_1, \mathcal{L}_2, \Theta)$, where $\text{op}$ is the operator name, $\mathcal{L}_1 = (\mathcal{G}_1, \mathcal{X}_1, \text{CRS}_1)$ is the primary vector layer with geometries $\mathcal{G}_1 = \{g_1, \dots, g_M\}$, attribute table $\mathcal{X}_1$, and coordinate reference system $\text{CRS}_1$, $\mathcal{L}_2$ is an optional secondary layer, and $\Theta$ is the dictionary of operator parameters.

### 3.1 Axiom 1: Metric and Geodesic CRS Invariant ($\mathcal{A}_{\text{CRS}}$)

Let $\mathcal{O}_{\text{metric}}$ be the set of operators that compute Euclidean distances, areas, spatial weights, or polygon overlays (`buffer`, `area`, `length`, `distance`, `nearest`, `areal_interpolation`, `spatial_join`, `dissolve`, `intersection`, `clip_and_area`, `morans_i`, `lisa`, `gwr`), and let $\mathcal{O}_{\text{area}} \subset \mathcal{O}_{\text{metric}}$ be the subset requiring equal-area preservation (`area`, `areal_interpolation`, `spatial_join`, `dissolve`). A workflow step $W$ satisfies $\mathcal{A}_{\text{CRS}}$ if and only if:

$$\text{op} \in \mathcal{O}_{\text{metric}} \implies \text{unit}(\text{CRS}_1) \in \{\text{m}, \text{ft}\} \land \neg\text{is\_geographic}(\text{CRS}_1), \qquad \mathcal{L}_2 \neq \varnothing \implies \text{CRS}_1 = \text{CRS}_2$$

**Distortion Calculation and Repair ($\mathcal{R}_{\text{CRS}}$):** When $\mathcal{L}_1$ is stored in `EPSG:4326` with bounding-box center $(\lambda_c, \phi_c)$ and spans $(\Delta\lambda, \Delta\phi)$, planar calculations in degree coordinates stretch East-West distances by $\sec(\phi_c) = 1/\cos(\phi_c)$. If $\mathcal{A}_{\text{CRS}}$ fails, $\mathcal{R}_{\text{CRS}}$ reprojects both $\mathcal{L}_1$ and $\mathcal{L}_2$ to $\text{CRS}^*$:
- `EPSG:5070` (NAD83 Conus Albers Equal-Area) if $(\text{op} \in \mathcal{O}_{\text{area}} \lor \Delta\lambda > 6.5^\circ)$ and the bounding box lies in CONUS;
- `EPSG:3035` (ETRS89 LAEA Europe) if $(\text{op} \in \mathcal{O}_{\text{area}} \lor \Delta\lambda > 6.5^\circ)$ and the bounding box lies in Europe;
- `EPSG:6933` (WGS 84 EASE-Grid 2.0 Global Equal-Area) for all other equal-area or wide-span operations;
- Local conformal UTM zone `EPSG:(32600 + Z)` (Northern Hemisphere) or `EPSG:(32700 + Z)` (Southern Hemisphere), where $Z = \lfloor (\lambda_c + 180)/6 \rfloor + 1 \in \{1,\dots,60\}$, for local distance, clip, and buffer operations.

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

### 3.3 Axiom 3: Topological 9-Intersection and Sliver Width Invariant ($\mathcal{A}_{\text{Topo}}$)

Every geometry $g \in \mathcal{G}$ must be non-empty and valid under the OGC Simple Features / Egenhofer 9-intersection rules (`g.is_valid == True`). For polygon overlay outputs (`intersection`, `difference`, `sym_difference`), we compare two sliver-pruning rules in projected metric coordinates:
1. **v1 Polsby-Popper + Relative Area Rule (`ns_vgp_v1`, `sliver_rule="pp_relative"`):** Prunes a polygon only when both $\text{PP}(g) = 4\pi A / P^2 < 0.005$ and $A / \overline{A}(\mathcal{G}) < 10^{-4}$. As shown in Section 4.4, coupling compactness with a relative area cap misses wider boundary slivers whose area exceeds $10^{-4}\overline{A}(\mathcal{G})$.
2. **v2 Mean-Width Rule (`ns_vgp`, `sliver_rule="mean_width"`):** For a narrow strip of length $L$ and width $w \ll L$, $A \approx L w$ and $P \approx 2(L+w) \approx 2L$, so the invariant mean width $W_{\text{mean}}(g) = 2A/P \approx w$ directly estimates the physical width of the strip in meters. Under `mean_width`, every output polygon in a projected metric CRS must satisfy:

$$W_{\text{mean}}(g) = \frac{2 \cdot \text{Area}(g)}{\text{Perimeter}(g)} \ge w_{\min} \quad (w_{\min} = 1.0\text{ m})$$

or satisfy the compactness and relative-area thresholds when $W_{\text{mean}}(g) \ge w_{\min}$.

### 3.4 Axiom 4: Spatial Statistical Pre-Conditions ($\mathcal{A}_{\text{Stat}}$)

For spatial statistical operators $\text{op} \in \{\text{morans\_i}, \text{lisa}, \text{getis\_ord}, \text{gwr}\}$ on target vector $\mathbf{y} \in \mathbb{R}^N$, spatial weights matrix $\mathbf{W} \in \mathbb{R}^{N \times N}$, and predictor matrix $\mathbf{X} \in \mathbb{R}^{N \times p}$:
1. $\mathbf{y}$ contains no `NaN` or $\pm\infty$ entries and $\text{Var}(\mathbf{y}) > 10^{-12}$.
2. Every spatial unit has at least one neighbor: $\forall i \in \{1,\dots,N\},\, \sum_{j=1}^N w_{ij} > 0$.
3. Sample size satisfies $N \ge 30$ for normal approximation, or conditional Monte Carlo permutation ($P \ge 999$) is used when $10 \le N < 30$.
4. For multivariate regression (`gwr`), the standardized design matrix $\tilde{\mathbf{X}}$ has condition number $\kappa(\tilde{\mathbf{X}}) < 30$.

---

## 4. Experimental Methodology, Independent Ground Truth, and Results

### 4.1 Scope of Baselines and Independent Ground-Truth Evaluator
To ensure full reproducibility without requiring paid API keys or non-deterministic LLM sampling, all baselines in `benchmarks/run_experiments.py` are deterministic GeoPandas reference implementations representing common code patterns:
- **`naive` (Unprojected GeoPandas Script):** Executes the spatial operator directly on the input `GeoDataFrame` using `1 deg = 111,320 m` for degree-based buffers and areas, `gpd.sjoin(..., predicate="intersects").groupby("zone_id").sum()` for spatial joins, `dissolve(by=..., aggfunc="sum")` for dissolves, and Queen contiguity weights (`weights="queen"`) without `NaN` imputation or collinearity pruning.
- **`naive_pre` (Standard Pre-Processing Baseline):** Applies standard geometry and table hygiene before running the operator: repairs invalid input geometries via `.make_valid()`, aligns secondary layer CRS to the primary layer (`secondary_gdf.to_crs(primary_gdf.crs)`), and drops `NaN` rows in statistical targets. It does not perform equal-area/UTM reprojection, semantic scale typing, sliver pruning, or weight-island healing.
- **`areaw` (Standard Area-Weighted Interpolation):** Projects both source and target polygons to `EPSG:5070` and interpolates extensive counts in proportion to exact overlap area $a_{s,t} / |A_s|$ and intensive variables by area-weighted averaging $a_{s,t} / \sum_{s'} a_{s',t}$, without Tobler's pycnophylactic smoothing.
- **`ns_vgp_v1` and `ns_vgp` (Ours):** Executes `NeuroSymbolicGeoVerifier(auto_heal=True)` with `sliver_rule="pp_relative"` (`ns_vgp_v1`) and `sliver_rule="mean_width"` (`ns_vgp`). Single-axiom ablations (`abl_no_A_CRS`, `abl_no_A_MAUP`, `abl_no_A_TOPO`, `abl_no_A_STAT`, `abl_all_off`) pass `disabled_axioms={...}` to `NeuroSymbolicGeoVerifier` so that disabling an axiom genuinely skips that axiom's verification and repair code path during execution.

Crucially, task pass/fail outcomes are **not** judged by NS-VGP's internal certificate flag. Instead, `benchmarks/run_experiments.py` evaluates every method's output `GeoDataFrame` against an **independent ground-truth oracle**:
1. **CRS Ground Truth (`CRS-01` to `CRS-15`):** True polygon areas are computed directly on the WGS84 ellipsoid (`a = 6378137.0 m`, `f = 1/298.257223563`) via `pyproj.Geod(ellps="WGS84").geometry_area_perimeter()`. True $1,000\text{ m}$ buffer and clipped areas are computed in a local Azimuthal Equidistant projection (`+proj=aeqd +lat_0=... +lon_0=... +datum=WGS84 +units=m`) centered on each task's centroid, which is distinct from both UTM and Albers projections used by NS-VGP. A run passes if and only if relative area error is $< 1.0\%$.
2. **MAUP Ground Truth (`MAUP-01` to `MAUP-15`):** For each task, we generate $40,000$ synthetic individual points across the study bounding box ($70\%$ drawn from $3$ bivariate Gaussian population hotspots and $30\%$ uniform background), each carrying an individual household count and local income. Aggregating the $40,000$ points into the $4 \times 4$ source tracts yields the observed input layer (`tracts`), while aggregating the exact same $40,000$ points into the $3 \times 3$ target zones (`zones`) yields the true target population and true mean income per zone. None of the polygon interpolation methods (`naive`, `naive_pre`, `areaw`, `ns_vgp`) ever see the $40,000$ points. A run passes if total population error is $< 0.5\%$, per-zone population NMAE is $< 15.0\%$, and intensive income MAPE is $< 10.0\%$.
3. **Topology Ground Truth (`TOPO-01` to `TOPO-15`):** Each task contains a self-intersecting bowtie parcel, two legitimate parcels, and $1$ to $3$ narrow boundary slivers with widths between $0.05\text{ m}$ and $0.60\text{ m}$ in `EPSG:32618`. The independent evaluator projects the output to `EPSG:32618` and computes the short side of each polygon's `minimum_rotated_rectangle`. Any output polygon with short-side width $< 1.0\text{ m}$ is counted as an unpruned sliver. A run passes if all geometries are valid, zero slivers remain, and $\ge 99.0\%$ of legitimate parcel area is retained.
4. **Spatial Statistics Ground Truth (`STAT-01` to `STAT-15`):** Each task has $N \in \{16, 20, 24\}$, one `NaN` in `disease_rate`, $2$ genuinely isolated island polygons shifted by $\pm 0.25^\circ$ (`dx = +0.25` for `k = N-1` and `dx = -0.25` for `k = N-2`), and a collinear predictor $x_2 = 2x_1 + \epsilon$ ($\kappa \in [61706, 133728]$). The evaluator inspects `out.attrs` to verify finite output, zero disconnected weight islands (`weights_islands == 0`), positive minimum row sum, and $\kappa(\tilde{\mathbf{X}}) < 30$.
5. **Clean Control Tasks (`CTRL-01` to `CTRL-10`):** Ten tasks with clean, already-projected `EPSG:32618` inputs (`buffer` and `area`) verify that NS-VGP does not corrupt valid workflows or trigger false-positive repairs.

All experiments were run on Windows with Python `3.11.9`, GeoPandas `1.1.4`, Shapely `2.2.0`, NumPy `1.26.4`, and Pandas `2.2.2` using `np.random.default_rng(42)`. Running `python benchmarks/run_experiments.py` executes all $535$ task-method combinations and writes `benchmarks/per_task_validation_log.csv` and `benchmarks/benchmark_results.json`.

![Figure 2: Benchmark Comparison](figures/fig2_benchmark_comparison.png)

| Task Category | Method | Passed ($k/N$) | Pass Rate (%) [95% Wilson CI] | Mean Area Error (%) | Total Pop Error (%) | Zone Pop NMAE (%) | Income Ratio ($\times$) | Mean Latency ($\text{ms} \pm \text{SD}$) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. CRS Distortion** ($N=15$) | `naive` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | $35.00\%$ (med $26.87\%$) | — | — | — | $58.00 \pm 166.33$ |
| | `naive_pre` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | $35.00\%$ (med $26.87\%$) | — | — | — | $14.04 \pm 8.05$ |
| | **`ns_vgp` (Ours)** | **$15/15$** | **$100.0\%$ $[79.6, 100.0]$** | **$0.023\%$ (med $0.015\%$)** | — | — | — | **$33.34 \pm 25.81$** |
| **2. MAUP Aggregation** ($N=15$) | `naive` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | — | $166.05\%$ (joins) | $166.05\%$ | $4.001\times$ | $23.76 \pm 13.48$ |
| | `naive_pre` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | — | $166.05\%$ (joins) | $166.05\%$ | $4.001\times$ | $19.44 \pm 8.51$ |
| | `areaw` ($N=10$ joins) | $10/10$ | $100.0\%$ $[72.2, 100.0]$ | — | $0.0003\%$ | $7.92\%$ | $1.001\times$ | $73.26 \pm 33.90$ |
| | **`ns_vgp` (Ours)** | **$15/15$** | **$100.0\%$ $[79.6, 100.0]$** | — | **$0.0000\%$ ($1.9\text{e-}15\%$)** | **$3.62\%$** | **$1.001\times$** | **$189.46 \pm 104.72$** |
| **3. Topology & Slivers** ($N=15$) | `naive` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | — | — | $2.00$ slivers | — | $235.47 \pm 31.87$ |
| | `naive_pre` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | — | — | $2.00$ slivers | — | $231.32 \pm 24.49$ |
| | `ns_vgp_v1` (`pp_relative`) | $10/15$ | $66.7\%$ $[41.7, 84.8]$ | — | — | $0.47$ slivers | — | $69.02 \pm 20.53$ |
| | **`ns_vgp` (`mean_width`)** | **$15/15$** | **$100.0\%$ $[79.6, 100.0]$** | — | — | **$0.00$ slivers** | — | **$66.13 \pm 11.89$** |
| **4. Spatial Statistics** ($N=15$) | `naive` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | — | — | $2$ islands, `NaN` | $\kappa \le 133,728$ | $7.71 \pm 1.48$ |
| | `naive_pre` | $0/15$ | $0.0\%$ $[0.0, 20.4]$ | — | — | $2$ islands | $\kappa \le 133,728$ | $10.65 \pm 1.98$ |
| | **`ns_vgp` (Ours)** | **$15/15$** | **$100.0\%$ $[79.6, 100.0]$** | — | — | **$0$ islands** | **$\kappa = 1.00$** | **$20.12 \pm 4.50$** |
| **5. Clean Control** ($N=10$) | `naive` | $10/10$ | $100.0\%$ $[72.2, 100.0]$ | $0.030\%$ | — | $0$ repairs | — | $5.65 \pm 5.04$ |
| | `naive_pre` | $10/10$ | $100.0\%$ $[72.2, 100.0]$ | $0.030\%$ | — | $0$ repairs | — | $5.65 \pm 5.36$ |
| | **`ns_vgp` (Ours)** | **$10/10$** | **$100.0\%$ $[72.2, 100.0]$** | **$0.030\%$** | — | **$0$ repairs** | — | **$13.55 \pm 3.85$** |
| **All Tasks ($N=70$)** | `naive` | $10/70$ | $14.3\%$ $[7.9, 24.3]$ | $35.00\%$ | $166.05\%$ | $166.05\%$ | $4.001\times$ | $70.56$ |
| | `naive_pre` | $10/70$ | $14.3\%$ $[7.9, 24.3]$ | $35.00\%$ | $166.05\%$ | $166.05\%$ | $4.001\times$ | $59.97$ |
| | **`ns_vgp` (Ours)** | **$70/70$** | **$100.0\%$ $[94.8, 100.0]$** | **$0.023\%$** | **$0.0000\%$** | **$3.62\%$** | **$1.001\times$** | **$68.16$** |

---

### 4.2 Granular Task-by-Task Validation Log (`per_task_validation_log.csv`)

The table below shows individual task records extracted directly from `benchmarks/per_task_validation_log.csv`:

| Task ID | Operator & Setup | Unverified Baseline (`naive` / `naive_pre`) | NS-VGP (`ns_vgp`) Output & Repair Trace | Error (`naive` $\to$ `ns_vgp`) | Latency (`ns_vgp`, ms) |
| :--- | :--- | :--- | :--- | :---: | :---: |
| `CRS-01` | `clip_and_area` ($\phi=12.0^\circ\text{N}$) | Estimate $=100,376,353\text{ m}^2$ vs Geod truth $=97,561,735\text{ m}^2$ | Reprojected layers to local UTM; estimate $=97,499,587\text{ m}^2$ ($1$ repair) | $2.88\% \to 0.064\%$ | $30.51$ |
| `CRS-06` | `area` ($\phi=37.8^\circ\text{N}$) | Estimate $=178,446,851\text{ m}^2$ vs Geod truth $=140,650,034\text{ m}^2$ | Reprojected to `EPSG:5070` Equal-Area; estimate $=140,650,031\text{ m}^2$ | $26.87\% \to 0.000002\%$ | $113.36$ |
| `CRS-07` | `clip_and_area` ($\phi=40.8^\circ\text{N}$) | Estimate $=100,376,353\text{ m}^2$ vs `aeqd` truth $=75,846,361\text{ m}^2$ | Reprojected both layers to `EPSG:32618`; estimate $=75,853,144\text{ m}^2$ | $32.34\% \to 0.009\%$ | $61.36$ |
| `CRS-12` | `area` ($\phi=55.7^\circ\text{N}$) | Estimate $=178,446,851\text{ m}^2$ vs Geod truth $=100,650,789\text{ m}^2$ | Reprojected to `EPSG:6933` Equal-Area; estimate $=100,650,794\text{ m}^2$ | $77.29\% \to 0.000004\%$ | $11.07$ |
| `CRS-14` | `area` ($\phi=60.2^\circ\text{N}$) | Estimate $=178,446,851\text{ m}^2$ vs Geod truth $=83,924,855\text{ m}^2$ | Reprojected to `EPSG:6933` Equal-Area; estimate $=83,924,860\text{ m}^2$ | $112.63\% \to 0.000005\%$ | $15.61$ |
| `CRS-15` | `buffer` ($1000\text{ m}$, $\phi=64.1^\circ\text{N}$) | Estimate $=160,615,171\text{ m}^2$ vs `aeqd` truth $=221,624,746\text{ m}^2$ | Reprojected to local UTM `EPSG:32628`; estimate $=221,668,055\text{ m}^2$ | $27.53\% \to 0.020\%$ | $46.46$ |
| `MAUP-01` | `spatial_join` ($16 \to 9$ zones) | True pop $=37,315$; `naive` $=105,618$ ($+183.04\%$); `areaw` zone NMAE $=10.63\%$ | Pycnophylactic pop $=37,315.0$ ($0.0\%$ total err, **$3.16\%$ zone NMAE**); income $\$84,938$ vs $\$84,900$ | $183.04\% \to 3.16\%$ NMAE | $305.68$ |
| `MAUP-04` | `spatial_join` ($16 \to 9$ zones) | `naive` pop err $=168.06\%$; `areaw` zone NMAE $=5.85\%$ | Pycnophylactic total err $=0.0\%$, **zone NMAE $=1.80\%$**; income MAPE $=0.77\%$ | $168.06\% \to 1.80\%$ NMAE | $214.17$ |
| `MAUP-03` | `dissolve` ($16 \to 4$ districts) | Summed intensive `median_income` ($4.002\times$ inflation, $300.30\%$ MAPE) | Blocked `SUM` on $\tau_{\text{int}}$; area-weighted income ratio $=1.0006\times$ ($1.15\%$ MAPE) | $300.30\% \to 1.15\%$ MAPE | $142.16$ |
| `MAUP-14` | `spatial_join` ($16 \to 9$ zones) | `naive` pop err $=163.21\%$; `areaw` zone NMAE $=7.49\%$ | Pycnophylactic total err $=0.0\%$, **zone NMAE $=1.64\%$**; income MAPE $=0.84\%$ | $163.21\% \to 1.64\%$ NMAE | $187.08$ |
| `TOPO-02` | `intersection` (bowtie + $2$ slivers) | `naive`/`naive_pre`: $2$ slivers; `ns_vgp_v1`: $1$ wider sliver unpruned | Repaired bowtie; `mean_width` ($2A/P < 1\text{ m}$) pruned all $2$ slivers; $100\%$ area kept | $2 \to 0$ slivers | $69.94$ |
| `STAT-01` | `gwr` ($N=16$, $2$ islands) | `naive`: `finite=False`, $2$ Queen islands, $\kappa=97,596.3$ | Imputed `NaN`; row-std $k$-NN ($0$ islands); pruned collinear column ($\kappa=1.00$) | Invalid $\to$ Valid | $17.40$ |
| `CTRL-01` | `buffer` (`EPSG:32618` clean) | Already projected; `naive` err $=0.038\%$ vs Geod | Zero repairs triggered (`repairs=0`); identical output to `naive` ($0.0\%$ diff) | $0.038\% \to 0.038\%$ | $20.59$ |

---

### 4.3 Cartographic Case Study (`MAUP-01`) and Point-Truth Comparison

![Figure 3: MAUP Case Study](figures/fig3_maup_case_study.png)

Figure 3 illustrates the point-truth evaluation on task `MAUP-01`:
- **Underlying Ground Truth:** $40,000$ synthetic individual residents with Gaussian hotspot clustering fall inside the study region ($37,315$ total population across the $16$ source tracts, with a true zone-level mean household income of $\$84,900$).
- **Naive Spatial Join (`naive` / `naive_pre`):** Running `gpd.sjoin(zones, tracts, predicate="intersects")` followed by `.groupby("zone_id").sum()` counts every boundary-straddling tract multiple times, inflating total population to $105,618$ ($+183.04\%$ error) and inflating mean zone income to $\$339,820$ ($4.00\times$ ecological fallacy).
- **Standard Area-Weighted Interpolation (`areaw`):** Projecting to `EPSG:5070` and weighting by overlap area $a_{s,t}/|A_s|$ conserves total population ($37,314.87$, $0.0004\%$ error), but because it assumes uniform population density inside each tract, its per-zone population NMAE against the true point counts is $10.63\%$ ($7.92\%$ mean across `MAUP-01..10`).
- **NS-VGP Pycnophylactic Interpolation (`ns_vgp`):** By reconstructing a smooth Laplacian density surface across tract boundaries while enforcing exact zone-wise mass rescaling, NS-VGP conserves total population exactly ($37,315.0$, $0.0000\%$ error) and reduces per-zone population NMAE on `MAUP-01` to **$3.16\%$** ($3.62\%$ mean across `MAUP-01..10`, beating uniform area weighting on $9$ of the $10$ spatial join tasks).

---

### 4.4 Ablation Analysis (`abl_no_A_*` and `ns_vgp_v1` vs `ns_vgp`)

![Figure 4: Ablation and Latency](figures/fig4_ablation_and_latency.png)

Figure 4(a) shows the empirical ablation study across the $60$ adversarial tasks (`CRS`, `MAUP`, `TOPO`, `STAT`), where each ablation disables a specific axiom inside `NeuroSymbolicGeoVerifier`:
- **Full `ns_vgp`:** Passes $60/60$ ($100.0\%$) of adversarial tasks and $70/70$ ($100.0\%$) overall.
- **No $\mathcal{A}_{\text{CRS}}$ (`abl_no_A_CRS`):** Passes $45/60$ ($75.0\%$), failing all $15$ `CRS` tasks ($0/15$) while passing `MAUP`, `TOPO`, and `STAT` ($15/15$ each) because their respective repair routines project layers locally when computing equal-area weights, metric sliver widths, and Euclidean $k$-NN centroids.
- **No $\mathcal{A}_{\text{MAUP}}$ (`abl_no_A_MAUP`):** Passes $45/60$ ($75.0\%$), failing all $15$ `MAUP` tasks ($0/15$) due to $166.05\%$ population duplication and $4.00\times$ income inflation.
- **No $\mathcal{A}_{\text{Topo}}$ (`abl_no_A_TOPO`):** Passes $45/60$ ($75.0\%$), failing all $15$ `TOPO` tasks ($0/15$) with $2.0$ unpruned slivers per task.
- **No $\mathcal{A}_{\text{Stat}}$ (`abl_no_A_STAT`):** Passes $45/60$ ($75.0\%$), failing all $15$ `STAT` tasks ($0/15$) due to un-imputed `NaN` values, $2$ Queen contiguity islands, and $\kappa \in [61706, 133728]$.
- **All Axioms Off (`abl_all_off`):** Passes $0/60$ ($0.0\%$).
- **Sliver Rule Ablation (`ns_vgp_v1` vs `ns_vgp`):** Under the independent minimum-rotated-rectangle evaluator (short-side width $< 1.0\text{ m}$), the v1 `pp_relative` rule passed $10/15$ (`TOPO` tasks, $66.7\%$) because 5 tasks contained wider boundary slivers ($0.35\text{--}0.60\text{ m}$ wide) whose area exceeded $10^{-4}\overline{A}(\mathcal{G})$. Replacing the relative-area cap with the invariant mean-width threshold $2A/P < 1.0\text{ m}$ in `ns_vgp` achieved $15/15$ ($100.0\%$) while retaining $100.0\%$ of legitimate parcel area.

Figure 4(b) shows the measured runtime distribution across categories: $13.55 \pm 3.85\text{ ms}$ on `CONTROL`, $20.12 \pm 4.50\text{ ms}$ on `STAT`, $33.34 \pm 25.81\text{ ms}$ on `CRS`, $66.13 \pm 11.89\text{ ms}$ on `TOPO`, and $189.46 \pm 104.72\text{ ms}$ on `MAUP` ($68.16\text{ ms}$ overall mean across all $70$ tasks, compared to $70.56\text{ ms}$ for `naive` and $59.97\text{ ms}$ for `naive_pre`). Notice that on `TOPO`, `ns_vgp` is faster than `naive` ($66.13\text{ ms}$ vs $235.47\text{ ms}$) because repairing the self-intersecting bowtie before running `gpd.overlay` avoids slow degenerate ring handling inside GEOS.

---

## 5. Limitations and Future Work

1. **Deterministic Reference Baselines vs. Live LLM Sampling:** Our benchmark evaluates deterministic GeoPandas implementations (`naive`, `naive_pre`, `areaw`) rather than live multi-turn LLM-Geo conversations. While these scripts capture the exact code patterns generated by LLMs, future evaluations should measure end-to-end closed-loop retry behavior across multiple commercial and open-weights LLMs.
2. **Lexical and Heuristic Column Typing:** The scale-type classifier $\tau(X_k)$ relies on regular-expression matching over column names combined with value-range checks. Ambiguously named columns (such as `val_1` or `attr_b`) without metadata units can still be misclassified unless explicit schema annotations or Core Concept Data Type (CCD) tags are provided in the layer metadata.
3. **Pycnophylactic Grid Resolution:** Our $40 \times 40$ Laplacian smoothing grid adds $\sim 189\text{ ms}$ per spatial join and assumes smooth density transitions across tract boundaries. Where sharp physical barriers (such as rivers, coastlines, or industrial zones) separate residential areas, dasymetric land-cover masks should be combined with pycnophylactic smoothing.

---

## 6. References

1. **Li, Z., & Ning, H. (2023).** Autonomous GIS: the next-generation AI-powered GIS. *International Journal of Digital Earth*, 16(2), 4668–4686. DOI: [10.1080/17538947.2023.2278895](https://doi.org/10.1080/17538947.2023.2278895)
2. **Akinboyewa, T., Li, Z., Ning, H., & Lessani, M. N. (2025).** GIS Copilot: towards an autonomous GIS agent for spatial analysis. *International Journal of Digital Earth*, 18(1). DOI: [10.1080/17538947.2025.2497489](https://doi.org/10.1080/17538947.2025.2497489)
3. **Ning, H., Li, Z., Akinboyewa, T., & Lessani, M. N. (2025).** An autonomous GIS agent framework for geospatial data retrieval. *International Journal of Digital Earth*, 18(1). DOI: [10.1080/17538947.2025.2458688](https://doi.org/10.1080/17538947.2025.2458688)
4. **Zhang, Y., Wei, C., He, Z., & Yu, W. (2024).** GeoGPT: An assistant for understanding and processing geospatial tasks. *International Journal of Applied Earth Observation and Geoinformation*, 131, 103976. DOI: [10.1016/j.jag.2024.103976](https://doi.org/10.1016/j.jag.2024.103976)
5. **Janowicz, K., Gao, S., McKenzie, G., Hu, Y., & Bhaduri, B. (2020).** GeoAI: spatially explicit artificial intelligence techniques for geographic knowledge discovery and beyond. *International Journal of Geographical Information Science*, 34(4), 625–636. DOI: [10.1080/13658816.2019.1684500](https://doi.org/10.1080/13658816.2019.1684500)
6. **Mai, G., Huang, W., Sun, J., Song, S., Mishra, D., et al. (2024).** On the Opportunities and Challenges of Foundation Models for GeoAI (Vision Paper). *ACM Transactions on Spatial Algorithms and Systems*, 10(2), Article 11, 1–46. DOI: [10.1145/3653070](https://doi.org/10.1145/3653070)
7. **Tobler, W. R. (1979).** Smooth Pycnophylactic Interpolation for Geographical Regions. *Journal of the American Statistical Association*, 74(367), 519–530. DOI: [10.1080/01621459.1979.10481647](https://doi.org/10.1080/01621459.1979.10481647)
8. **Openshaw, S. (1984).** *The Modifiable Areal Unit Problem*. CATMOG No. 38, Geo Books.
9. **Goodchild, M. F., & Lam, N. S.-N. (1980).** Areal interpolation: a variant of the traditional spatial problem. *Geo-Processing*, 1(3), 297–312.
10. **Kuhn, W. (2012).** Core concepts of spatial information for transdisciplinary research. *International Journal of Geographical Information Science*, 26(12), 2267–2276. DOI: [10.1080/13658816.2012.722637](https://doi.org/10.1080/13658816.2012.722637)
11. **Scheider, S., & Ballatore, A. (2017).** Semantic typing of linked geoprocessing workflows. *International Journal of Digital Earth*, 11(1), 113–138. DOI: [10.1080/17538947.2017.1305457](https://doi.org/10.1080/17538947.2017.1305457)
12. **Scheider, S., Meerlo, R., Kasalica, V., & Lamprecht, A.-L. (2020).** Ontology of core concept data types for answering geo-analytical questions. *Journal of Spatial Information Science*, 20, 167–201. DOI: [10.5311/JOSIS.2020.20.620](https://doi.org/10.5311/JOSIS.2020.20.620)
13. **Egenhofer, M. J., & Franzosa, R. D. (1991).** Point-set topological spatial relations. *International Journal of Geographical Information Systems*, 5(2), 161–174. DOI: [10.1080/02693799108927841](https://doi.org/10.1080/02693799108927841)
14. **Polsby, D. D., & Popper, R. D. (1991).** The Third Criterion: Compactness as a Procedural Safeguard Against Partisan Gerrymandering. *Yale Law & Policy Review*, 9(2), 301–353.
15. **Anselin, L. (1995).** Local Indicators of Spatial Association—LISA. *Geographical Analysis*, 27(2), 93–115. DOI: [10.1111/j.1538-4632.1995.tb00338.x](https://doi.org/10.1111/j.1538-4632.1995.tb00338.x)
16. **Brunsdon, C., Fotheringham, A. S., & Charlton, M. E. (1996).** Geographically Weighted Regression: A Method for Exploring Spatial Nonstationarity. *Geographical Analysis*, 28(4), 281–298. DOI: [10.1111/j.1538-4632.1996.tb00936.x](https://doi.org/10.1111/j.1538-4632.1996.tb00936.x)
17. **d'Avila Garcez, A., & Lamb, L. C. (2023).** Neurosymbolic AI: The 3rd Wave. *Artificial Intelligence Review*, 56, 12387–12406. DOI: [10.1007/s10462-023-10448-w](https://doi.org/10.1007/s10462-023-10448-w)
