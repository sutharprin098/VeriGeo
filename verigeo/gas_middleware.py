"""GAS 1.0.0 & Multi-Step DAG Verification Middleware for VeriGeo (NS-VGP 1.0).

Provides:
  1. `VeriGeoGasMiddleware`: Intercepts Geospatial Agentic Services (GAS 1.0.0)
     `ExecuteTask` payloads, verifies and auto-heals spatial pre/post-conditions,
     and attaches a cryptographic `SpatialProofCertificate` (`spc_certificate.json`)
     to the GAS task response `artifacts` and `diagnostics`.
  2. `VerifiedDAGExecutor`: Executes multi-step spatial analysis Directed Acyclic
     Graphs (DAGs) where every node transition is verified across the 4 formal
     GIScience axioms (A_CRS, A_MAUP, A_TOPO, A_STAT), producing a chained
     merkle-style workflow proof certificate.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

import geopandas as gpd

from .protocol import (
    PROTOCOL_VERSION,
    NeuroSymbolicGeoVerifier,
    SpatialProofCertificate,
)


class VeriGeoGasMiddleware:
    """Neuro-Symbolic Verification Middleware for GAS 1.0.0 Agent Services."""

    def __init__(self, auto_heal: bool = True) -> None:
        self.verifier = NeuroSymbolicGeoVerifier(auto_heal=auto_heal)

    def execute_verified_gas_task(
        self,
        agent_id: str,
        operator: str,
        primary_gdf: gpd.GeoDataFrame,
        secondary_gdf: gpd.GeoDataFrame | None = None,
        parameters: dict[str, Any] | None = None,
        semantic_overrides: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """Execute a GAS 1.0.0 task with formal 4-axiom verification and SPC attachment."""
        result_gdf, spc = self.verifier.verify_and_execute(
            operator=operator,
            primary_gdf=primary_gdf,
            secondary_gdf=secondary_gdf,
            parameters=parameters,
            semantic_overrides=semantic_overrides,
        )
        spc_dict = spc.to_dict()
        return {
            "service": "GAS",
            "version": "1.0.0",
            "verigeo_protocol": PROTOCOL_VERSION,
            "agent_id": agent_id,
            "status": "SUCCEEDED" if spc.verifiable_guarantee else "FAILED_VERIFICATION",
            "verifiable_guarantee": spc.verifiable_guarantee,
            "scientific_validity_score": spc.scientific_validity_score,
            "result_gdf": result_gdf,
            "spatial_proof_certificate": spc_dict,
            "diagnostics": {
                "repairs_applied": spc.repairs_applied,
                "provenance_sha256": spc.provenance_sha256,
                "execution_latency_ms": spc.execution_latency_ms,
            },
        }


class VerifiedDAGExecutor:
    """Multi-Step Spatial Workflow DAG Executor with Chained Spatial Proof Certificates."""

    def __init__(self, auto_heal: bool = True) -> None:
        self.verifier = NeuroSymbolicGeoVerifier(auto_heal=auto_heal)

    def execute_dag(
        self,
        nodes: list[dict[str, Any]],
        initial_layers: dict[str, gpd.GeoDataFrame],
    ) -> dict[str, Any]:
        """Execute a topological sequence of spatial workflow nodes with formal verification.

        Each node dict specifies:
          - `node_id`: Unique step identifier
          - `operator`: Spatial operator name
          - `primary_input`: Key in `initial_layers` or previous `node_id`
          - `secondary_input`: Optional key in `initial_layers` or previous `node_id`
          - `parameters`: Optional dict of operator parameters
        """
        state: dict[str, Any] = dict(initial_layers)
        certificates: list[SpatialProofCertificate] = []
        chained_hashes: list[str] = []

        for node in nodes:
            node_id = str(node["node_id"])
            op = str(node["operator"])
            p_key = str(node["primary_input"])
            s_key = node.get("secondary_input")
            params = dict(node.get("parameters") or {})

            gdf1 = state[p_key]
            gdf2 = state[s_key] if s_key and s_key in state else None

            out_gdf, spc = self.verifier.verify_and_execute(
                operator=op,
                primary_gdf=gdf1,
                secondary_gdf=gdf2,
                parameters=params,
            )
            state[node_id] = out_gdf
            certificates.append(spc)
            chained_hashes.append(spc.provenance_sha256)

        dag_merkle_root = hashlib.sha256(
            "|".join(chained_hashes).encode("utf-8")
        ).hexdigest()

        return {
            "protocol_version": PROTOCOL_VERSION,
            "num_nodes": len(nodes),
            "all_nodes_verified": all(c.verifiable_guarantee for c in certificates),
            "mean_scientific_validity_score": (
                sum(c.scientific_validity_score for c in certificates) / max(len(certificates), 1)
            ),
            "dag_merkle_root_sha256": dag_merkle_root,
            "node_certificates": [c.to_dict() for c in certificates],
            "outputs": state,
        }
