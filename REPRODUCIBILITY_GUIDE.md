# VeriGeo Independent Reproducibility Guide

This guide enables researchers to reproduce all benchmark executions, figures, tables, and verification records from source code without proprietary API keys.

## 1. Environment Setup

Ensure Python 3.10+ is installed. Clone the repository and install dependencies in editable mode:

```bash
git clone https://github.com/sutharprin098/VeriGeo.git
cd VeriGeo
pip install -e .[dev]
```

Dependencies:
- `geopandas >= 0.14.0`
- `shapely >= 2.0.0`
- `pyproj >= 3.6.0`
- `numpy >= 1.24.0`
- `pandas >= 2.0.0`
- `matplotlib >= 3.7.0`
- `pytest >= 7.0.0`

## 2. Running the Test Suite

Execute the 9 unit and integration tests:

```bash
pytest -v
```

Previously recorded output: `8 passed in ~1.5s` on the original suite. The revision branch adds `test_stat13_iterative_condition_precheck_only()`; its result is pending and must be verified with a fresh run.

## 3. Running the Complete 70-Task Benchmark Suite (535 Executions)

To execute all $70$ benchmark tasks across the 3 main methods (`naive`, `naive_pre`, `ns_vgp`), the area-weighted baseline (`areaw`), the v1 compactness rule (`ns_vgp_v1`), and the five ablation configurations ($5 \times 60$ runs):

```bash
python benchmarks/run_experiments.py
```

Outputs generated:
- `benchmarks/per_task_validation_log.csv`: Per-task log with 535 rows recording execution metrics, errors, and outcomes.
- `benchmarks/benchmark_results.json`: Summary statistics, pass rates with 95% Wilson confidence intervals, and aggregate error metrics.
- `paper/figures/fig1_protocol_architecture.png` / `.pdf`: Pipeline schematic.
- `paper/figures/fig2_benchmark_comparison.png` / `.pdf`: Pass rates and log-scale error reduction.
- `paper/figures/fig3_maup_case_study.png` / `.pdf`: 40,000-point ground truth vs. spatial join methods.
- `paper/figures/fig4_ablation_and_latency.png` / `.pdf`: Ablation matrix and per-category latency breakdown.

## 4. Rebuilding the LaTeX Manuscript and Overleaf Bundle

To verify LaTeX manuscript compilation and regenerate the Overleaf submission zip:

```bash
python -c "
import zipfile, pathlib
p = pathlib.Path('paper')
zf = zipfile.ZipFile(p / 'VeriGeo_Overleaf_Submission.zip', 'w', zipfile.ZIP_DEFLATED)
zf.write(p / 'main.tex', 'main.tex')
zf.write(p / 'references.bib', 'references.bib')
for f in sorted((p / 'figures').glob('*')):
    zf.write(f, f'figures/{f.name}')
zf.close()
print('Overleaf bundle rebuilt successfully.')
"
```
