# GlueTK Validation & Benchmark

This folder contains simple, reproducible benchmarks used in the JCIM manuscript.

## Quick start (PyMOL)

From the PyMOL command line:

```python
run validation/benchmark_analysis.py
benchmark_all()
```

Outputs:
- `validation/results/benchmark_summary.csv`
- Per-case CSVs under `validation/results/cases/`

## Customize dataset

Edit `BENCHMARK_CASES` in `validation/benchmark_analysis.py` to match the PDB IDs, chain IDs, and (optionally) glue resnames used in your manuscript.

