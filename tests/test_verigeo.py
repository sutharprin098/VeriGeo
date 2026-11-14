"""Unit and integration tests for VeriGeo (NS-VGP 1.0): Neuro-Symbolic Verifiable GeoAI Protocol."""

from __future__ import annotations

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, box

from verigeo import (
    PROTOCOL_VERSION,
    NeuroSymbolicGeoVerifier,
    VeriGeoGasMiddleware,
    VerifiedDAGExecutor,
    classify_attribute_semantics,
    compute_optimal_metric_crs,
    polsby_popper_compactness,
    pycnophylactic_areal_interpolation,
    tool_verify_and_heal_spatial_workflow,
)


def test_classify_attribute_semantics_ontology() -> None:
    """Verify semantic type classification distinguishes extensive, intensive, and nominal variables."""
    assert classify_attribute_semantics("population") == "extensive"
    assert classify_attribute_semantics("total_households") == "extensive"
    assert classify_attribute_semantics("emissions_tons") == "extensive"
    assert classify_attribute_semantics("pop_density") == "intensive"
    assert classify_attribute_semantics("median_income") == "intensive"
    assert classify_attribute_semantics("crime_rate_per_1k") == "intensive"
    assert classify_attribute_semantics("ndvi") == "intensive"
    assert classify_attribute_semantics("pm25") == "intensive"
    assert classify_attribute_semantics("tract_id") == "nominal"
    assert classify_attribute_semantics("geoid") == "nominal"


def test_axiom1_crs_metric_reprojection_and_distortion() -> None:
    """Verify Axiom 1 detects unprojected geographic CRS and auto-reprojects to metric UTM/Equal-Area."""
    gdf = gpd.GeoDataFrame(
        {"id": ["A", "B"], "population": [1000, 2000]},
        geometry=[box(-77.88, 40.78, -77.84, 40.82), box(-77.84, 40.78, -77.80, 40.82)],
        crs="EPSG:4326",
    )
    verifier = NeuroSymbolicGeoVerifier(auto_heal=True)
    res_gdf, spc = verifier.verify_and_execute(
        operator="buffer",
        primary_gdf=gdf,
        parameters={"distance_m": 500.0},
    )
    assert spc.verifiable_guarantee is True
    assert spc.scientific_validity_score == 1.0
    crs_ax = next(a for a in spc.axiom_results if a.axiom_id == "A_CRS")
    assert crs_ax.status == "AUTO_REPAIRED"
    assert crs_ax.pre_condition_met is False
    assert crs_ax.post_condition_met is True
    assert not res_gdf.crs.is_geographic
    assert str(res_gdf.crs).startswith("EPSG:326")


def test_axiom2_maup_pycnophylactic_mass_conservation_and_intensive_scale() -> None:
    """Verify Axiom 2 blocks SUM on intensive variables and conserves exact mass on extensive counts."""
    src_polys = [
        box(-77.90, 40.75, -77.85, 40.80),
        box(-77.85, 40.75, -77.80, 40.80),
        box(-77.90, 40.80, -77.85, 40.85),
        box(-77.85, 40.80, -77.80, 40.85),
    ]
    src_gdf = gpd.GeoDataFrame(
        {
            "tract_id": ["T1", "T2", "T3", "T4"],
            "population": [2500.0, 4500.0, 3000.0, 5000.0],
            "median_income": [50000.0, 70000.0, 60000.0, 80000.0],
            "ndvi": [0.35, 0.55, 0.45, 0.65],
        },
        geometry=src_polys,
        crs="EPSG:4326",
    )
    tgt_polys = [
        box(-77.90, 40.75, -77.83, 40.82),
        box(-77.83, 40.75, -77.80, 40.82),
        box(-77.90, 40.82, -77.80, 40.85),
    ]
    tgt_gdf = gpd.GeoDataFrame(
        {"zone_id": ["Z1", "Z2", "Z3"]},
        geometry=tgt_polys,
        crs="EPSG:4326",
    )

    verifier = NeuroSymbolicGeoVerifier(auto_heal=True)
    res_gdf, spc = verifier.verify_and_execute(
        operator="spatial_join",
        primary_gdf=src_gdf,
        secondary_gdf=tgt_gdf,
        parameters={"aggfunc": "sum"},
    )
    assert spc.verifiable_guarantee is True
    maup_ax = next(a for a in spc.axiom_results if a.axiom_id == "A_MAUP")
    assert maup_ax.status == "AUTO_REPAIRED"
    assert maup_ax.post_condition_met is True

    # Extensive mass must be 100% conserved (15,000 total population)
    assert abs(float(res_gdf["population"].sum()) - 15000.0) < 1e-3
    # Intensive median_income must stay within [50000, 80000] (never summed to >150,000!)
    assert float(res_gdf["median_income"].min()) >= 50000.0 - 1.0
    assert float(res_gdf["median_income"].max()) <= 80000.0 + 1.0


def test_axiom3_topology_self_intersection_and_sliver_pruning() -> None:
    """Verify Axiom 3 repairs bowtie self-intersections and prunes Polsby-Popper overlay slivers."""
    bowtie = Polygon(
        [
            (-77.85, 40.78),
            (-77.81, 40.82),
            (-77.81, 40.78),
            (-77.85, 40.82),
            (-77.85, 40.78),
        ]
    )
    valid_poly = box(-77.81, 40.78, -77.77, 40.82)
    sliver_poly = box(-77.77, 40.78, -77.7699995, 40.82)

    gdf1 = gpd.GeoDataFrame(
        {"id": ["bowtie", "valid", "sliver"]},
        geometry=[bowtie, valid_poly, sliver_poly],
        crs="EPSG:4326",
    )
    gdf2 = gpd.GeoDataFrame(
        {"study": ["S1"]},
        geometry=[box(-77.86, 40.77, -77.76, 40.83)],
        crs="EPSG:4326",
    )

    verifier = NeuroSymbolicGeoVerifier(auto_heal=True)
    res_gdf, spc = verifier.verify_and_execute(
        operator="intersection",
        primary_gdf=gdf1,
        secondary_gdf=gdf2,
    )
    assert spc.verifiable_guarantee is True
    topo_ax = next(a for a in spc.axiom_results if a.axiom_id == "A_TOPO")
    assert topo_ax.status == "AUTO_REPAIRED"
    assert topo_ax.post_condition_met is True
    assert all(res_gdf.geometry.is_valid)
    assert len(res_gdf) == 2


def test_axiom4_spatial_statistics_island_and_collinearity_healing() -> None:
    """Verify Axiom 4 heals NaN values, disconnected Queen islands, small N, and GWR multicollinearity."""
    polys = [box(-77.90 + i * 0.03, 40.75, -77.87 + i * 0.03, 40.78) for i in range(10)]
    polys.append(box(-77.10, 41.20, -77.07, 41.23))
    n = len(polys)
    x1 = np.linspace(10.0, 25.0, n)
    x2_collinear = 3.0 * x1 + 1e-5
    y = 2.0 * x1 + 5.0
    y[1] = np.nan

    gdf = gpd.GeoDataFrame(
        {"rate": y, "x1": x1, "x2": x2_collinear},
        geometry=polys,
        crs="EPSG:4326",
    )
    verifier = NeuroSymbolicGeoVerifier(auto_heal=True)
    res_gdf, spc = verifier.verify_and_execute(
        operator="gwr",
        primary_gdf=gdf,
        parameters={"column": "rate", "weights": "queen", "predictors": ["x1", "x2"]},
    )
    assert spc.verifiable_guarantee is True
    stat_ax = next(a for a in spc.axiom_results if a.axiom_id == "A_STAT")
    assert stat_ax.status == "AUTO_REPAIRED"
    assert stat_ax.post_condition_met is True
    assert stat_ax.metrics["contiguity_islands_detected"] >= 1
    assert stat_ax.metrics["final_weights_scheme"] == "knn_row_standardized"
    assert "global_moran_i" in res_gdf.columns


def test_gas_middleware_and_verified_dag_executor() -> None:
    """Verify GAS 1.0.0 middleware and multi-node VerifiedDAGExecutor produce Merkle SPC proofs."""
    src_gdf = gpd.GeoDataFrame(
        {
            "tract_id": ["T1", "T2"],
            "population": [3000.0, 7000.0],
            "median_income": [55000.0, 85000.0],
        },
        geometry=[box(-77.90, 40.75, -77.85, 40.80), box(-77.85, 40.75, -77.80, 40.80)],
        crs="EPSG:4326",
    )
    tgt_gdf = gpd.GeoDataFrame(
        {"zone_id": ["Z1", "Z2"]},
        geometry=[box(-77.90, 40.75, -77.84, 40.80), box(-77.84, 40.75, -77.80, 40.80)],
        crs="EPSG:4326",
    )

    mw = VeriGeoGasMiddleware(auto_heal=True)
    gas_resp = mw.execute_verified_gas_task(
        agent_id="vector_analysis_agent",
        operator="spatial_join",
        primary_gdf=src_gdf,
        secondary_gdf=tgt_gdf,
        parameters={"aggfunc": "sum"},
    )
    assert gas_resp["status"] == "SUCCEEDED"
    assert gas_resp["verifiable_guarantee"] is True

    dag = VerifiedDAGExecutor(auto_heal=True)
    dag_out = dag.execute_dag(
        nodes=[
            {
                "node_id": "step1_interp",
                "operator": "spatial_join",
                "primary_input": "tracts",
                "secondary_input": "zones",
                "parameters": {"aggfunc": "sum"},
            },
            {
                "node_id": "step2_buffer",
                "operator": "buffer",
                "primary_input": "step1_interp",
                "parameters": {"distance_m": 250.0},
            },
        ],
        initial_layers={"tracts": src_gdf, "zones": tgt_gdf},
    )
    assert dag_out["all_nodes_verified"] is True
    assert dag_out["mean_scientific_validity_score"] == 1.0
    assert len(dag_out["dag_merkle_root_sha256"]) == 64


def test_tool_handler_verify_and_heal() -> None:
    """Verify the tool wrapper returns a valid Spatial Proof Certificate."""
    out = tool_verify_and_heal_spatial_workflow(
        iface=None,
        args={
            "operator": "spatial_join",
            "parameters": {"aggfunc": "sum"},
            "auto_heal": True,
        },
    )
    assert out["ok"] is True
    assert out["protocol"] == PROTOCOL_VERSION
    assert out["verifiable_guarantee"] is True
    assert out["scientific_validity_score"] == 1.0
