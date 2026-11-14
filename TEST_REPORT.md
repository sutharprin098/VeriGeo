# Independent Test and Validation Report

## 1. Automated Test Suite Execution
- **Command**: `pytest -v d:\VeriGeo\tests`
- **Environment**: Python 3.11.9, GeoPandas 1.1.4, Shapely 2.2.0, NumPy 1.26.4, Pandas 2.2.2
- **Platform**: Windows 11 (AMD64)
- **Status**: PASSED (8 passed, 0 failed in 1.29s)

### Test Case Results:
1. `test_classify_attribute_semantics_ontology`: PASSED
   - Verified that `population`, `total_households`, `emissions_tons` are typed as `extensive`.
   - Verified that `pop_density`, `median_income`, `crime_rate_per_1k`, `ndvi`, `pm25` are typed as `intensive`.
   - Verified that `tract_id`, `geoid` are typed as `nominal`.
2. `test_axiom1_crs_metric_reprojection_and_distortion`: PASSED
   - Geographic CRS (`EPSG:4326`) buffer input automatically reprojected to local UTM (`EPSG:32618`).
   - Verified that output CRS is metric and distortion is eliminated.
3. `test_axiom2_maup_pycnophylactic_mass_conservation_and_intensive_scale`: PASSED
   - Tested $16 \to 3$ polygon aggregation.
   - Total extensive population conserved to exact value ($15,000.0$).
   - Intensive `median_income` preserved within $[50000, 80000]$ source range (never summed to $>150,000$).
4. `test_axiom3_topology_self_intersection_and_sliver_pruning`: PASSED
   - Bowtie self-intersection repaired via `make_valid`.
   - Overlay sliver with width $< 1.0\text{ m}$ pruned.
   - Legitimate parcel retained.
5. `test_axiom4_spatial_statistics_island_and_collinearity_healing`: PASSED
   - `NaN` imputed with spatial median.
   - Queen contiguity island healed via metric $k$-NN fallback.
   - Multicollinear predictor pruned with ridge regularization enabled.
6. `test_service_middleware_and_verified_dag_executor`: PASSED
   - Multi-node DAG (spatial join $\to$ buffer) executed with Merkle root hash verification.
7. `test_tool_handler_verify_and_heal`: PASSED
   - Standalone tool interface emits valid Spatial Verification Record (SPC).
8. `test_disabled_axioms_ablation_and_sliver_rules`: PASSED
   - Verified that `disabled_axioms={'A_CRS'}` skips projection repair for ablation benchmarking.

## 2. Sensitivity Analysis Results: Pycnophylactic Grid Resolution vs. Iterations
Evaluated on $4 \to 2$ incongruent polygon transfer (analytical uniform area-weighted benchmark: $Z_1 = 7000.0, Z_2 = 3000.0$):
- **$20 \times 20$ grid**:
  - $15$ iters: $Z_1 = 6835.06, Z_2 = 3164.94$, Sum $= 10000.0000$
  - $35$ iters: $Z_1 = 6779.45, Z_2 = 3220.55$, Sum $= 10000.0000$
  - $70$ iters: $Z_1 = 6754.03, Z_2 = 3245.97$, Sum $= 10000.0000$
- **$40 \times 40$ grid (Default)**:
  - $15$ iters: $Z_1 = 6983.16, Z_2 = 3016.84$, Sum $= 10000.0000$
  - $35$ iters: $Z_1 = 6935.86, Z_2 = 3064.14$, Sum $= 10000.0000$
  - $70$ iters: $Z_1 = 6881.21, Z_2 = 3118.79$, Sum $= 10000.0000$
- **$80 \times 80$ grid**:
  - $35$ iters: $Z_1 = 6997.24, Z_2 = 3002.76$, Sum $= 10000.0000$
- **$160 \times 160$ grid**:
  - $35$ iters: $Z_1 = 7000.00, Z_2 = 3000.00$, Sum $= 10000.0000$

**Finding**: Mass is 100.0000% conserved across all grid resolutions and iteration counts. Finer grids converge toward exact analytical polygon boundaries with smooth intra-zone density modulation.
