"""Measured benchmark for VeriGeo (NS-VGP 1.0).

Everything reported here is computed by running code. Nothing is hard-coded.

Methods compared
----------------
naive      : idiomatic one-liner GeoPandas scripts, i.e. what an LLM or a thin service wrapper
             typically emits (1 deg = 111,320 m, sjoin + groupby.sum, plain overlay, queen weights).
             NOTE: these are our own scripts. They are NOT live external LLM agents.
naive_pre  : naive + light preprocessing (CRS alignment, make_valid, dropna).
areaw      : plain area-weighted interpolation without pycnophylactic smoothing (MAUP tasks only).
ns_vgp     : VeriGeo with all four axioms enabled.
abl_no_*   : VeriGeo with one axiom truly disabled. abl_all_off: all four disabled.

Ground truth is independent of VeriGeo's own checks
---------------------------------------------------
CRS    : geodesic area on the WGS84 ellipsoid (pyproj.Geod) and a local azimuthal-equidistant
         projection for buffers and clips.
MAUP   : a synthetic population of individual points; source tracts and target zones are both
         aggregated from the SAME points, so the target values are exact.
TOPO   : validity from shapely, and slivers defined by mean polygon width (2A/P) < 1 m,
         which differs from VeriGeo's own Polsby-Popper rule.
STAT   : finiteness of outputs, zero-row-sum count of the weights actually used, and the
         condition number of the predictors actually used.

Run:  python benchmarks/run_experiments.py
"""

from __future__ import annotations

import csv
import json
import math
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from pyproj import CRS, Geod
from shapely.geometry import Polygon, box
from shapely.validation import make_valid

warnings.filterwarnings("ignore")

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from verigeo import NeuroSymbolicGeoVerifier  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent
FIG_DIR = ROOT_DIR / "paper" / "figures"
GEOD = Geod(ellps="WGS84")
KM = 111320.0  # metres per degree assumed by naive scripts
ALL_AXIOMS = frozenset({"A_CRS", "A_MAUP", "A_TOPO", "A_STAT"})
MAUP_COLS = ["population", "households", "pop_density", "median_income", "pm25", "ndvi"]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def geo_area(geom: Any) -> float:
    if geom is None or geom.is_empty:
        return 0.0
    return abs(GEOD.geometry_area_perimeter(geom)[0])


def geo_area_sum(gdf: gpd.GeoDataFrame) -> float:
    return float(sum(geo_area(g) for g in gdf.to_crs("EPSG:4326").geometry))


def aeqd(lon: float, lat: float) -> CRS:
    return CRS.from_proj4(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs")


def to_metric(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    if gdf.crs is not None and gdf.crs.is_geographic:
        return gdf.to_crs(gdf.estimate_utm_crs())
    return gdf


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    r = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - r) / d) * 100, min(1.0, (c + r) / d) * 100)


def rel_err(est: float, truth: float) -> float:
    if not np.isfinite(est) or truth == 0:
        return float("inf")
    return abs(est - truth) / abs(truth) * 100.0


# ----------------------------------------------------------------------------
# Task construction
# ----------------------------------------------------------------------------
LATS = [12.0, 18.5, 25.0, 28.6, 34.0, 37.8, 40.8, 43.7, 45.5, 48.9, 51.5, 55.7, 59.9, 62.0, 64.1]


def build_crs_tasks() -> list[dict[str, Any]]:
    rng = np.random.default_rng(42)
    tasks = []
    for i in range(15):
        lat = LATS[i]
        lon = -77.8 + (i * 7.5) % 120.0 - 40.0
        polys = [
            box(lon + dx * 0.04, lat + dy * 0.04, lon + (dx + 1) * 0.04, lat + (dy + 1) * 0.04)
            for dx in range(3)
            for dy in range(3)
        ]
        g1 = gpd.GeoDataFrame(
            {"feature_id": [f"F_{k}" for k in range(9)], "population": rng.integers(1200, 8500, 9)},
            geometry=polys,
            crs="EPSG:4326",
        )
        op = "buffer" if i % 2 == 0 else "area"
        g2 = None
        if i % 3 == 0:
            op = "clip_and_area"
            g2 = gpd.GeoDataFrame(
                {"clip_id": ["C1"]},
                geometry=[box(lon + 0.01, lat + 0.01, lon + 0.10, lat + 0.10)],
                crs="EPSG:4326",
            ).to_crs("EPSG:3857")
        tasks.append(_crs_task(f"CRS-{i+1:02d}", "CRS", op, g1, g2, lat, lon))
    return tasks


def build_control_tasks() -> list[dict[str, Any]]:
    """Clean, already-projected inputs. A verifier must NOT change the answer here."""
    tasks = []
    for i in range(10):
        lat, lon = 40.8, -77.8 + 0.2 * i
        utm = CRS.from_epsg(32618)
        c = gpd.GeoSeries([shapely.geometry.Point(lon, lat)], crs=4326).to_crs(utm).iloc[0]
        polys = [
            box(c.x + dx * 800, c.y + dy * 800, c.x + (dx + 1) * 800, c.y + (dy + 1) * 800)
            for dx in range(3)
            for dy in range(3)
        ]
        g1 = gpd.GeoDataFrame({"feature_id": [f"F_{k}" for k in range(9)]}, geometry=polys, crs=utm)
        op = "buffer" if i % 2 == 0 else "area"
        tasks.append(_crs_task(f"CTRL-{i+1:02d}", "CONTROL", op, g1, None, lat, lon))
    return tasks


def _crs_task(tid, cat, op, g1, g2, lat, lon) -> dict[str, Any]:
    g4326 = g1.to_crs("EPSG:4326")
    if op == "area":
        truth = geo_area_sum(g4326)
    elif op == "buffer":
        total = 0.0
        for geom in g4326.geometry:
            ctr = geom.centroid
            local = gpd.GeoSeries([geom], crs=4326).to_crs(aeqd(ctr.x, ctr.y))
            buf = local.buffer(1000.0)
            total += geo_area_sum(gpd.GeoDataFrame(geometry=buf.to_crs(4326), crs=4326))
        truth = total
    else:
        loc = aeqd(lon + 0.055, lat + 0.055)
        a = g1.to_crs(loc)
        b = g2.to_crs(loc)
        truth = float(gpd.clip(a, b).geometry.area.sum())
    return {
        "task_id": tid, "category": cat, "operator": op, "gdf1": g1, "gdf2": g2,
        "params": {"distance_m": 1000.0, "buffer_distance": 1000.0},
        "meta": {"lat": lat, "lon": lon, "truth": truth},
    }


def build_maup_tasks() -> list[dict[str, Any]]:
    tasks = []
    x0, y0, W = -77.90, 40.75, 0.10
    for i in range(15):
        rng = np.random.default_rng(1000 + i)
        n = 40000
        centers = rng.uniform([x0 + 0.2 * W, y0 + 0.2 * W], [x0 + 0.8 * W, y0 + 0.8 * W], size=(3, 2))
        sig = rng.uniform(0.012, 0.03, 3)
        comp = rng.choice(4, size=n, p=[0.25, 0.25, 0.25, 0.25])  # 3 hotspots + uniform background
        pts = np.where(
            (comp < 3)[:, None],
            centers[np.minimum(comp, 2)] + rng.normal(size=(n, 2)) * sig[np.minimum(comp, 2), None],
            rng.uniform([x0, y0], [x0 + W, y0 + W], size=(n, 2)),
        )
        ok = (pts[:, 0] >= x0) & (pts[:, 0] < x0 + W) & (pts[:, 1] >= y0) & (pts[:, 1] < y0 + W)
        pts = pts[ok]
        px, py = pts[:, 0], pts[:, 1]
        fx, fy = (px - x0) / W, (py - y0) / W
        income = 45000 + 80000 * fx + rng.normal(0, 8000, len(px))
        pm25 = 8 + 12 * fy + rng.normal(0, 1.5, len(px))
        ndvi = np.clip(0.2 + 0.5 * (1 - fx) + rng.normal(0, 0.05, len(px)), 0, 1)

        def agg(ix, iy, nx, ny):
            cid = iy * nx + ix
            cnt = np.bincount(cid, minlength=nx * ny).astype(float)
            out = {"pop": cnt}
            for k, v in (("inc", income), ("pm", pm25), ("nd", ndvi)):
                s = np.bincount(cid, weights=v, minlength=nx * ny)
                out[k] = np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)
            return out

        ix4 = np.clip((fx * 4).astype(int), 0, 3)
        iy4 = np.clip((fy * 4).astype(int), 0, 3)
        T = agg(ix4, iy4, 4, 4)
        src_polys = [box(x0 + c * 0.025, y0 + r * 0.025, x0 + (c + 1) * 0.025, y0 + (r + 1) * 0.025)
                     for r in range(4) for c in range(4)]
        area_km2 = np.array([geo_area(p) for p in src_polys]) / 1e6
        src = gpd.GeoDataFrame(
            {
                "tract_id": [f"T_{k+1:02d}" for k in range(16)],
                "district": [f"D_{k // 4}" for k in range(16)],
                "population": T["pop"],
                "households": np.round(T["pop"] * 0.4),
                "pop_density": T["pop"] / area_km2,
                "median_income": T["inc"],
                "pm25": T["pm"],
                "ndvi": T["nd"],
            },
            geometry=src_polys,
            crs="EPSG:4326",
        )
        src = src.fillna({c: float(np.nanmean(src[c])) for c in ["median_income", "pm25", "ndvi"]})

        dissolve = i % 3 == 2
        if not dissolve:
            ix3 = np.clip((fx * 3).astype(int), 0, 2)
            iy3 = np.clip((fy * 3).astype(int), 0, 2)
            Z = agg(ix3, iy3, 3, 3)
            step = W / 3
            tgt = gpd.GeoDataFrame(
                {"zone_id": [f"Z_{k+1:02d}" for k in range(9)]},
                geometry=[box(x0 + c * step, y0 + r * step, x0 + (c + 1) * step, y0 + (r + 1) * step)
                          for r in range(3) for c in range(3)],
                crs="EPSG:4326",
            )
            truth = {"pop": Z["pop"], "inc": Z["inc"]}
            op, params = "spatial_join", {"aggfunc": "sum"}
        else:
            D = agg(np.zeros_like(iy4), iy4, 1, 4)
            tgt = None
            truth = {"pop": D["pop"], "inc": D["inc"]}
            op, params = "dissolve", {"aggfunc": "sum", "by": "district"}
        tasks.append({
            "task_id": f"MAUP-{i+1:02d}", "category": "MAUP", "operator": op,
            "gdf1": src, "gdf2": tgt, "params": params,
            "meta": {"truth": truth, "n_points": int(len(px))},
        })
    return tasks


def build_topo_tasks() -> list[dict[str, Any]]:
    tasks = []
    for i in range(15):
        rng = np.random.default_rng(2000 + i)
        bx = -77.85 + rng.uniform(-0.02, 0.02)
        by = 40.78 + rng.uniform(-0.02, 0.02)
        s = rng.uniform(0.03, 0.05)
        bow = Polygon([(bx, by), (bx + s, by + s), (bx + s, by), (bx, by + s), (bx, by)])
        p1 = box(bx + s, by, bx + 2 * s, by + s)
        p2 = box(bx, by + s, bx + 2 * s, by + 2 * s)
        slivers = []
        for k in range(1 + i % 3):
            w = rng.uniform(1e-7, 1.5e-6)
            x = bx + 2 * s + 0.004 * (k + 1)
            slivers.append(box(x, by, x + w, by + 2 * s))
        geoms = [bow, p1, p2] + slivers
        g1 = gpd.GeoDataFrame(
            {"parcel_id": [f"P{k}" for k in range(len(geoms))], "pop_total": np.ones(len(geoms))},
            geometry=geoms, crs="EPSG:4326",
        )
        g2 = gpd.GeoDataFrame(
            {"study_area": ["SA1"]},
            geometry=[box(bx - 0.005, by - 0.005, bx + 2 * s + 0.03, by + 2 * s + 0.005)],
            crs="EPSG:4326",
        )
        legit = gpd.GeoDataFrame(geometry=[make_valid(bow), p1, p2], crs=4326)
        legit_area = float(to_metric(legit).geometry.area.sum())
        tasks.append({
            "task_id": f"TOPO-{i+1:02d}", "category": "TOPO", "operator": "intersection",
            "gdf1": g1, "gdf2": g2, "params": {},
            "meta": {"legit_area": legit_area, "n_slivers_in": len(slivers)},
        })
    return tasks


def build_stat_tasks() -> list[dict[str, Any]]:
    tasks = []
    for i in range(15):
        rng = np.random.default_rng(3000 + i)
        n = 16 + (i % 3) * 4
        polys = []
        for k in range(n):
            r, c = k // 4, k % 4
            # two genuinely isolated single-polygon islands (zero neighbours under queen contiguity)
            dx = 0.25 if k == n - 1 else (-0.25 if k == n - 2 else 0.0)
            dy = 0.25 if k >= n - 2 else 0.0
            x0 = -77.90 + c * 0.03 + dx
            y0 = 40.75 + r * 0.03 + dy
            polys.append(box(x0, y0, x0 + 0.03, y0 + 0.03))
        vals = 50 + 4 * (np.arange(n) % 4) + rng.normal(0, 5, n)
        vals[int(rng.integers(0, n - 2))] = np.nan
        x1 = rng.normal(10, 2.5, n)
        x2 = 2.0 * x1 + rng.normal(0, 1e-4, n)
        op = ["gwr", "lisa", "morans_i"][i % 3]
        g = gpd.GeoDataFrame({"disease_rate": vals, "x1": x1, "x2_collinear": x2},
                             geometry=polys, crs="EPSG:4326")
        params: dict[str, Any] = {"column": "disease_rate", "weights": "queen"}
        if op == "gwr":
            params["predictors"] = ["x1", "x2_collinear"]
        tasks.append({"task_id": f"STAT-{i+1:02d}", "category": "STAT", "operator": op,
                      "gdf1": g, "gdf2": None, "params": params, "meta": {"n": n}})
    return tasks


# ----------------------------------------------------------------------------
# Runners + independent scoring (one function per category)
# ----------------------------------------------------------------------------
def _ns(method, op, g1, g2, params):
    rule = method[2] if len(method) > 2 else "mean_width"
    v = NeuroSymbolicGeoVerifier(auto_heal=True, disabled_axioms=method[1], sliver_rule=rule)
    return v.verify_and_execute(op, g1, g2, params)


def score_crs(task, method) -> dict[str, Any]:
    g1, g2, op, truth = task["gdf1"], task["gdf2"], task["operator"], task["meta"]["truth"]
    t0 = time.perf_counter()
    repairs = 0
    try:
        if method[0] in ("naive", "naive_pre"):
            geographic = g1.crs.is_geographic
            if op == "area":
                est = float(g1.geometry.area.sum() * (KM**2 if geographic else 1.0))
            elif op == "buffer":
                d = 1000.0 / KM if geographic else 1000.0
                est = geo_area_sum(gpd.GeoDataFrame(geometry=g1.geometry.buffer(d), crs=g1.crs))
            else:
                clipped = gpd.clip(g1, g2.to_crs(g1.crs))
                est = float(clipped.geometry.area.sum() * (KM**2 if geographic else 1.0))
        else:
            res, spc = _ns(method, op, g1, g2, dict(task["params"]))
            repairs = len(spc.repairs_applied)
            est = geo_area_sum(res) if op == "buffer" else float(res["area_m2"].sum())
        err = rel_err(est, truth)
        crashed = False
    except Exception as exc:  # noqa: BLE001
        est, err, crashed = float("nan"), float("inf"), True
    return {
        "passed": bool(err < 1.0), "crashed": crashed, "rel_area_err_pct": err,
        "estimate": est, "truth": truth, "repairs": repairs,
        "latency_ms": (time.perf_counter() - t0) * 1000,
    }


def _maup_arrays(task, method):
    src, tgt, op = task["gdf1"], task["gdf2"], task["operator"]
    if method[0] in ("naive", "naive_pre"):
        if op == "spatial_join":
            j = gpd.sjoin(tgt, src, how="left", predicate="intersects")
            a = j.groupby(j.index)[MAUP_COLS].sum().reindex(tgt.index)
        else:
            a = src[MAUP_COLS + ["district", "geometry"]].dissolve(by="district", aggfunc="sum")
            a = a.sort_index()
        return a["population"].to_numpy(float), a["median_income"].to_numpy(float), None
    if method[0] == "areaw":
        s, t = src.to_crs(5070).reset_index(drop=True), tgt.to_crs(5070).reset_index(drop=True)
        s["_sid"], t["_tid"] = np.arange(len(s)), np.arange(len(t))
        ov = gpd.overlay(s[["_sid"] + MAUP_COLS + ["geometry"]], t[["_tid", "geometry"]], how="intersection")
        ov["a"] = ov.geometry.area
        sa = s.geometry.area.to_numpy()
        ov["frac"] = ov["a"] / sa[ov["_sid"].to_numpy()]
        pop = (ov["population"] * ov["frac"]).groupby(ov["_tid"]).sum().reindex(range(len(t))).fillna(0)
        inc = ((ov["median_income"] * ov["a"]).groupby(ov["_tid"]).sum()
               / ov["a"].groupby(ov["_tid"]).sum()).reindex(range(len(t)))
        return pop.to_numpy(float), inc.to_numpy(float), None
    res, spc = _ns(method, task["operator"], src, tgt, dict(task["params"]))
    key = "zone_id" if task["operator"] == "spatial_join" else "district"
    r = res.sort_values(key)
    return r["population"].to_numpy(float), r["median_income"].to_numpy(float), spc


def score_maup(task, method) -> dict[str, Any] | None:
    if method[0] == "areaw" and task["operator"] != "spatial_join":
        return None
    t0 = time.perf_counter()
    truth = task["meta"]["truth"]
    try:
        pop, inc, spc = _maup_arrays(task, method)
        tp, ti = truth["pop"], truth["inc"]
        pop_total_err = rel_err(float(np.nansum(pop)), float(tp.sum()))
        pop_zone_nmae = float(np.nansum(np.abs(pop - tp)) / tp.sum() * 100)
        inc_ratio = float(np.nanmean(inc) / np.nanmean(ti))
        inc_mape = float(np.nanmean(np.abs(inc - ti) / ti) * 100)
        passed = bool(pop_total_err < 0.1 and 0.9 <= inc_ratio <= 1.1)
        crashed = False
    except Exception:  # noqa: BLE001
        pop = inc = np.array([])
        pop_total_err = pop_zone_nmae = inc_mape = float("inf")
        inc_ratio = float("nan")
        passed, crashed = False, True
    return {
        "passed": passed, "crashed": crashed, "pop_total_err_pct": pop_total_err,
        "pop_zone_nmae_pct": pop_zone_nmae, "income_mean_ratio": inc_ratio,
        "income_mape_pct": inc_mape, "latency_ms": (time.perf_counter() - t0) * 1000,
        "_pop": pop, "_inc": inc,
    }


def _short_side(g) -> float:
    """Short side (m) of the minimum rotated bounding rectangle. Independent of VeriGeo's rules."""
    try:
        rect = g.minimum_rotated_rectangle
        xs, ys = rect.exterior.coords.xy
        s1 = math.hypot(xs[1] - xs[0], ys[1] - ys[0])
        s2 = math.hypot(xs[2] - xs[1], ys[2] - ys[1])
        return min(s1, s2)
    except Exception:  # noqa: BLE001  (degenerate line/point)
        return 0.0


def topo_eval(out: gpd.GeoDataFrame, legit_area: float) -> dict[str, Any]:
    invalid = int((~out.geometry.is_valid).sum())
    parts = to_metric(out).explode(index_parts=False)
    parts = parts[~parts.geometry.is_empty]
    thin = parts.geometry.apply(_short_side) < 1.0
    return {
        "invalid_geoms": invalid,
        "slivers": int(thin.sum()),
        "legit_area_retained_pct": float(parts.area[~thin].sum() / legit_area * 100),
    }


def score_topo(task, method) -> dict[str, Any]:
    g1, g2 = task["gdf1"], task["gdf2"]
    t0 = time.perf_counter()
    crashed, out = False, None
    try:
        if method[0] == "naive":
            out = gpd.overlay(g1, g2, how="intersection", keep_geom_type=False)
        elif method[0] == "naive_pre":
            a = g1.copy()
            a["geometry"] = a.geometry.apply(lambda g: make_valid(g) if not g.is_valid else g)
            out = gpd.overlay(a, g2, how="intersection", keep_geom_type=False)
        else:
            out, _ = _ns(method, "intersection", g1, g2, {})
    except Exception:  # noqa: BLE001
        crashed = True
    if out is None:
        m = {"invalid_geoms": -1, "slivers": -1, "legit_area_retained_pct": 0.0}
        passed = False
    else:
        m = topo_eval(out, task["meta"]["legit_area"])
        passed = bool(m["invalid_geoms"] == 0 and m["slivers"] == 0)
    return {"passed": passed, "crashed": crashed, **m, "latency_ms": (time.perf_counter() - t0) * 1000}


def score_stat(task, method) -> dict[str, Any]:
    g, op = task["gdf1"], task["operator"]
    params = dict(task["params"])
    t0 = time.perf_counter()
    crashed = False
    try:
        if method[0] == "naive":
            v = NeuroSymbolicGeoVerifier(auto_heal=False, disabled_axioms=ALL_AXIOMS)
            res, _ = v.verify_and_execute(op, g, None, params)
        elif method[0] == "naive_pre":
            v = NeuroSymbolicGeoVerifier(auto_heal=False, disabled_axioms=ALL_AXIOMS)
            res, _ = v.verify_and_execute(op, g.dropna(subset=["disease_rate"]), None, params)
        else:
            res, _ = _ns(method, op, g, None, params)
        finite = bool(np.all(np.isfinite(res["local_moran_i"])) and np.isfinite(res["global_moran_i"].iloc[0]))
        islands = int(res.attrs.get("weights_islands", -1))
        kappa = 1.0
        preds = [p for p in res.attrs.get("predictors_used", []) if p in g.columns]
        if op == "gwr" and len(preds) >= 2:
            X = g[preds].to_numpy(float)
            X = (X - X.mean(0)) / np.clip(X.std(0), 1e-12, None)
            kappa = float(np.linalg.cond(X))
        passed = bool(finite and islands == 0 and kappa < 30.0)
    except Exception:  # noqa: BLE001
        finite, islands, kappa, passed, crashed = False, -1, float("nan"), False, True
    return {"passed": passed, "crashed": crashed, "finite_output": finite, "weight_islands": islands,
            "design_condition_number": kappa, "latency_ms": (time.perf_counter() - t0) * 1000}


SCORERS = {"CRS": score_crs, "CONTROL": score_crs, "MAUP": score_maup, "TOPO": score_topo, "STAT": score_stat}

METHODS = {
    "naive": ("naive",),
    "naive_pre": ("naive_pre",),
    "areaw": ("areaw",),
    "ns_vgp": ("ns", frozenset()),
    "ns_vgp_v1": ("ns", frozenset(), "pp_relative"),
    "abl_no_A_CRS": ("ns", frozenset({"A_CRS"})),
    "abl_no_A_MAUP": ("ns", frozenset({"A_MAUP"})),
    "abl_no_A_TOPO": ("ns", frozenset({"A_TOPO"})),
    "abl_no_A_STAT": ("ns", frozenset({"A_STAT"})),
    "abl_all_off": ("ns", ALL_AXIOMS),
}


# ----------------------------------------------------------------------------
# Figures (all values come from `rows` / `summary`)
# ----------------------------------------------------------------------------
def fig_architecture() -> None:
    fig, ax = plt.subplots(figsize=(12, 5.4))
    ax.set_xlim(0, 12), ax.set_ylim(0, 5.4), ax.axis("off")

    def b(x, y, w, h, t, s, fc, ec):
        ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1,rounding_size=0.15",
                                             fc=fc, ec=ec, lw=1.6))
        ax.text(x + w / 2, y + h * 0.68, t, ha="center", va="center", fontsize=9.5, fontweight="bold")
        ax.text(x + w / 2, y + h * 0.30, s, ha="center", va="center", fontsize=7.8)

    b(0.3, 3.3, 2.5, 1.4, "Planner (LLM / Agent Client)", "emits operator + parameters", "#eff6ff", "#2563eb")
    b(0.3, 0.6, 2.5, 1.6, "Input layers", "any CRS, any topology,\nmixed attribute types", "#f8fafc", "#64748b")
    ax.add_patch(mpatches.FancyBboxPatch((3.3, 0.4), 5.4, 4.6, boxstyle="round,pad=0.1,rounding_size=0.2",
                                         fc="#fefce8", ec="#ca8a04", lw=1.8, ls="--"))
    ax.text(6.0, 4.75, "Symbolic verifier + repair", ha="center", fontsize=10.5, fontweight="bold", color="#854d0e")
    rows = [
        ("A_CRS: metric CRS", "reproject to UTM / equal-area", 3.55),
        ("A_MAUP: scale types", "no SUM on intensive; mass-preserving transfer", 2.6),
        ("A_TOPO: valid geometry", "make_valid; drop slivers after overlay", 1.65),
        ("A_STAT: pre-conditions", "NaN, islands, small N, collinearity", 0.7),
    ]
    for t, s, y in rows:
        b(3.55, y, 4.9, 0.8, t, s, "#ffffff", "#d97706")
    b(9.2, 3.2, 2.5, 1.5, "Execute operator", "on repaired inputs", "#ecfdf5", "#059669")
    b(9.2, 0.6, 2.5, 1.7, "Spatial Proof\nCertificate", "pre/post flags, repairs,\nSHA-256 of the trace", "#f0fdf4", "#15803d")
    kw = dict(arrowstyle="->", lw=1.6, color="#334155")
    ax.annotate("", xy=(3.3, 4.0), xytext=(2.8, 4.0), arrowprops=kw)
    ax.annotate("", xy=(3.3, 1.4), xytext=(2.8, 1.4), arrowprops=kw)
    ax.annotate("", xy=(9.2, 3.9), xytext=(8.7, 3.9), arrowprops=kw)
    ax.annotate("", xy=(10.45, 2.3), xytext=(10.45, 3.2), arrowprops=dict(arrowstyle="<-", lw=1.6, color="#059669"))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig1_protocol_architecture.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_main(summary) -> None:
    cats = ["CRS", "MAUP", "TOPO", "STAT", "CONTROL"]
    names = {"naive": "Naive", "naive_pre": "Naive + preprocessing", "ns_vgp": "VeriGeo (NS-VGP)"}
    cols = {"naive": "#f87171", "naive_pre": "#fbbf24", "ns_vgp": "#10b981"}
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    x = np.arange(len(cats) + 1)
    w = 0.26
    for j, m in enumerate(["naive", "naive_pre", "ns_vgp"]):
        vals, lo, hi = [], [], []
        for c in cats + ["ALL"]:
            s = summary["pass_rates"][c][m]
            vals.append(s["pass_pct"]), lo.append(max(0.0, s["pass_pct"] - s["ci95"][0])), hi.append(max(0.0, s["ci95"][1] - s["pass_pct"]))
        axes[0].bar(x + (j - 1) * w, vals, w, yerr=[lo, hi], capsize=2.5, color=cols[m], label=names[m], ec="#334155")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(cats + ["ALL"])
    axes[0].set_ylabel("Tasks passing independent check (%)")
    axes[0].set_title("(a) Pass rate by category (Wilson 95% CI)")
    axes[0].set_ylim(0, 112)
    axes[0].grid(axis="y", ls=":", alpha=0.6)
    axes[0].legend(loc="upper right", fontsize=8)

    e = summary["error_metrics"]
    labels = ["CRS area\nerror (%)", "MAUP total\npop. error (%)", "MAUP income\nmean error (%)"]
    keys = ["crs_area_err_median", "maup_pop_total_err_mean", "maup_income_err_mean"]
    for j, m in enumerate(["naive", "ns_vgp"]):
        vals = [max(e[m][k], 1e-6) for k in keys]
        bars = axes[1].bar(np.arange(3) + (j - 0.5) * 0.36, vals, 0.36, color=cols[m], label=names[m], ec="#334155")
        for bar, v in zip(bars, [e[m][k] for k in keys]):
            axes[1].text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.15, f"{v:.3g}",
                         ha="center", fontsize=8)
    axes[1].set_yscale("log")
    axes[1].set_xticks(range(3))
    axes[1].set_xticklabels(labels)
    axes[1].set_ylabel("Measured error (%) (log scale)")
    axes[1].set_title("(b) Measured error vs independent ground truth")
    axes[1].grid(axis="y", ls=":", alpha=0.6)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig2_benchmark_comparison.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_case(tasks) -> dict[str, Any]:
    t = [x for x in tasks if x["task_id"] == "MAUP-01"][0]
    truth = t["meta"]["truth"]
    npop, ninc, _ = _maup_arrays(t, ("naive",))
    spop, sinc, _ = _maup_arrays(t, ("ns", frozenset()))
    apop, ainc, _ = _maup_arrays(t, ("areaw",))
    tgt = t["gdf2"].copy()
    vmax = max(truth["pop"].max(), npop.max(), spop.max())
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9))
    for ax, (title, vals) in zip(axes.ravel()[:3], [
        (f"(a) Ground truth (from points)\ntotal = {truth['pop'].sum():,.0f}", truth["pop"]),
        (f"(b) Naive sjoin + sum\ntotal = {npop.sum():,.0f}", npop),
        (f"(c) VeriGeo\ntotal = {spop.sum():,.1f}", spop),
    ]):
        tgt["v"] = vals
        tgt.plot(column="v", ax=ax, cmap="YlGnBu", edgecolor="#1e293b", vmin=0, vmax=vmax, legend=True)
        ax.set_xticks([-77.89, -77.85, -77.81])
        ax.set_title(title, fontsize=10)
    ax = axes[1, 1]
    ax.scatter(truth["pop"], npop, c="#f87171", label="Naive", s=36, ec="k")
    ax.scatter(truth["pop"], apop, c="#fbbf24", label="Area-weighted", s=36, ec="k", marker="s")
    ax.scatter(truth["pop"], spop, c="#10b981", label="VeriGeo", s=36, ec="k", marker="^")
    lim = [0, vmax * 1.05]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xlabel("True zone population"), ax.set_ylabel("Estimated zone population")
    ax.set_title("(d) Per-zone estimate vs truth")
    ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig3_maup_case_study.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {
        "task": "MAUP-01", "true_total_pop": float(truth["pop"].sum()),
        "naive_total_pop": float(npop.sum()), "verigeo_total_pop": float(spop.sum()),
        "areaw_total_pop": float(apop.sum()),
        "true_mean_income": float(np.nanmean(truth["inc"])), "naive_mean_income": float(np.nanmean(ninc)),
        "verigeo_mean_income": float(np.nanmean(sinc)),
        "zone_nmae_pct": {
            "naive": float(np.abs(npop - truth["pop"]).sum() / truth["pop"].sum() * 100),
            "areaw": float(np.abs(apop - truth["pop"]).sum() / truth["pop"].sum() * 100),
            "verigeo": float(np.abs(spop - truth["pop"]).sum() / truth["pop"].sum() * 100),
        },
    }


def fig_ablation(summary) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5))
    order = ["ns_vgp", "abl_no_A_CRS", "abl_no_A_MAUP", "abl_no_A_TOPO", "abl_no_A_STAT", "abl_all_off"]
    lab = ["Full", "no A_CRS", "no A_MAUP", "no A_TOPO", "no A_STAT", "all off"]
    cats = ["CRS", "MAUP", "TOPO", "STAT"]
    mat = np.array([[summary["pass_rates"][c][m]["pass_pct"] for c in cats] for m in order])
    im = axes[0].imshow(mat, vmin=0, vmax=100, cmap="RdYlGn", aspect="auto")
    axes[0].set_xticks(range(4)), axes[0].set_xticklabels(cats)
    axes[0].set_yticks(range(len(order))), axes[0].set_yticklabels(lab)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            axes[0].text(j, i, f"{mat[i, j]:.0f}", ha="center", va="center", fontsize=9)
    axes[0].set_title("(a) Real ablation: pass rate (%) with each axiom disabled")
    fig.colorbar(im, ax=axes[0], fraction=0.04)
    lat = summary["latency_ms"]
    names = ["CRS", "MAUP", "TOPO", "STAT", "CONTROL"]
    axes[1].bar(names, [lat[n]["mean"] for n in names], yerr=[lat[n]["std"] for n in names],
                capsize=4, color="#3b82f6", ec="#1d4ed8")
    for i, n in enumerate(names):
        axes[1].text(i, lat[n]["mean"] + lat[n]["std"] + 3, f"{lat[n]['mean']:.1f}", ha="center", fontsize=8)
    axes[1].set_ylabel("VeriGeo wall-clock per task (ms)")
    axes[1].set_title("(b) Measured runtime (mean +/- SD)")
    axes[1].grid(axis="y", ls=":", alpha=0.6)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig4_ablation_and_latency.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> dict[str, Any]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    tasks = (build_crs_tasks() + build_maup_tasks() + build_topo_tasks()
             + build_stat_tasks() + build_control_tasks())
    rows: list[dict[str, Any]] = []
    for t in tasks:
        for mname, spec in METHODS.items():
            if mname.startswith("abl_") and t["category"] == "CONTROL":
                continue
            if mname == "areaw" and t["category"] != "MAUP":
                continue
            if mname == "ns_vgp_v1" and t["category"] != "TOPO":
                continue
            r = SCORERS[t["category"]](t, spec)
            if r is None:
                continue
            clean = {k: v for k, v in r.items() if not k.startswith("_")}
            rows.append({"task_id": t["task_id"], "category": t["category"],
                         "operator": t["operator"], "method": mname, **clean})
        print(t["task_id"], "done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "per_task_validation_log.csv", index=False)

    cats = ["CRS", "MAUP", "TOPO", "STAT", "CONTROL"]
    pass_rates: dict[str, Any] = {}
    for c in cats + ["ALL"]:
        sub = df if c == "ALL" else df[df.category == c]
        pass_rates[c] = {}
        for m in df.method.unique():
            s = sub[sub.method == m]
            if len(s) == 0:
                continue
            k, n = int(s.passed.sum()), int(len(s))
            pass_rates[c][m] = {"k": k, "n": n, "pass_pct": k / n * 100, "ci95": list(wilson(k, n))}

    def agg(cat, m, col, how="mean", where=None):
        s = df[(df.category == cat) & (df.method == m)]
        if where is not None:
            s = s[s.operator == where]
        v = s[col].replace([np.inf, -np.inf], np.nan).dropna()
        return float(getattr(v, how)()) if len(v) else float("nan")

    error_metrics = {}
    for m in ["naive", "ns_vgp"]:
        s = df[(df.category == "MAUP") & (df.method == m) & (df.operator == "spatial_join")]
        error_metrics[m] = {
            "crs_area_err_median": agg("CRS", m, "rel_area_err_pct", "median"),
            "crs_area_err_mean": agg("CRS", m, "rel_area_err_pct", "mean"),
            "maup_pop_total_err_mean": float(s.pop_total_err_pct.mean()),
            "maup_pop_zone_nmae_mean": float(s.pop_zone_nmae_pct.mean()),
            "maup_income_ratio_mean": float(s.income_mean_ratio.mean()),
            "maup_income_err_mean": float((s.income_mean_ratio - 1).abs().mean() * 100),
            "maup_income_mape_mean": float(s.income_mape_pct.mean()),
        }
    sj = df[(df.category == "MAUP") & (df.method == "areaw")]
    error_metrics["areaw"] = {
        "maup_pop_total_err_mean": float(sj.pop_total_err_pct.mean()),
        "maup_pop_zone_nmae_mean": float(sj.pop_zone_nmae_pct.mean()),
        "maup_income_mape_mean": float(sj.income_mape_pct.mean()),
    }
    ctrl = df[df.category == "CONTROL"]
    piv = ctrl.pivot(index="task_id", columns="method", values="estimate")
    control_check = {
        "max_rel_diff_ns_vs_naive": float(((piv["ns_vgp"] - piv["naive"]).abs() / piv["naive"]).max()),
        "ns_repairs_total": int(ctrl[ctrl.method == "ns_vgp"].repairs.sum()),
    }
    topo = {m: {
        "crashes": int(df[(df.category == "TOPO") & (df.method == m)].crashed.sum()),
        "mean_slivers": float(df[(df.category == "TOPO") & (df.method == m)].slivers.clip(lower=0).mean()),
        "mean_invalid": float(df[(df.category == "TOPO") & (df.method == m)].invalid_geoms.clip(lower=0).mean()),
        "mean_legit_area_retained_pct": float(
            df[(df.category == "TOPO") & (df.method == m)].legit_area_retained_pct.mean()),
    } for m in ["naive", "naive_pre", "ns_vgp_v1", "ns_vgp"]}
    latency = {}
    for c in cats:
        s = df[(df.category == c) & (df.method == "ns_vgp")].latency_ms
        latency[c] = {"mean": float(s.mean()), "std": float(s.std(ddof=0))}
    base_lat = {m: float(df[(df.method == m)].latency_ms.mean()) for m in ["naive", "naive_pre", "ns_vgp"]}

    summary = {
        "environment": {"python": sys.version.split()[0], "geopandas": gpd.__version__,
                        "shapely": shapely.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        "n_tasks": {c: int(sum(1 for t in tasks if t["category"] == c)) for c in cats},
        "n_rows_logged": int(len(df)),
        "pass_rates": pass_rates, "error_metrics": error_metrics, "control_check": control_check,
        "topology": topo, "latency_ms": latency, "mean_latency_by_method_ms": base_lat,
    }
    fig_architecture()
    fig_main(summary)
    summary["case_study"] = fig_case(tasks)
    fig_ablation(summary)
    (OUT_DIR / "benchmark_results.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, default=float))
