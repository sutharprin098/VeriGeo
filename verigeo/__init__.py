"""VeriGeo (NS-VGP 1.0): Neuro-Symbolic Verifiable GeoAI Protocol."""

from .protocol import (
    EQUAL_AREA_OPERATORS,
    METRIC_REQUIRED_OPERATORS,
    PROTOCOL_VERSION,
    SPATIAL_STAT_OPERATORS,
    AxiomCheckResult,
    NeuroSymbolicGeoVerifier,
    SpatialProofCertificate,
    classify_attribute_semantics,
    compute_optimal_metric_crs,
    polsby_popper_compactness,
    pycnophylactic_areal_interpolation,
    tool_verify_and_heal_spatial_workflow,
)
from .middleware import VeriGeoServiceMiddleware, VerifiedDAGExecutor

__version__ = "1.0.0"

__all__ = [
    "__version__",
    "PROTOCOL_VERSION",
    "METRIC_REQUIRED_OPERATORS",
    "EQUAL_AREA_OPERATORS",
    "SPATIAL_STAT_OPERATORS",
    "AxiomCheckResult",
    "SpatialProofCertificate",
    "NeuroSymbolicGeoVerifier",
    "VeriGeoServiceMiddleware",
    "VerifiedDAGExecutor",
    "classify_attribute_semantics",
    "compute_optimal_metric_crs",
    "polsby_popper_compactness",
    "pycnophylactic_areal_interpolation",
    "tool_verify_and_heal_spatial_workflow",
]
