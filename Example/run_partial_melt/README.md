# Partial-Melt Workflow

This directory contains the staged partial-melt entry point. It starts from a
full-melt GCE case, freezes the core, and steps the silicate reservoir through a
chain of lower melt fractions.

For the standard full-melt pipeline, use [`../README.md`](../README.md).

## Python Environment

From the repository root (`GlobalChemicalEquilibrium/`), install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Run commands from the repository root with `PYTHONPATH=.` so imports such as
`tools` and `Example` resolve correctly.

## Run

```bash
PYTHONPATH=. python3 Example/run_partial_melt/run_partial_melt.py
```

Before running, edit [`run_partial_melt.py`](run_partial_melt.py). Most runs
only need:

- `params`: the `PartialMeltParams(...)` science setup.
- `run_name`: the label used for the partial-melt results directory.
- `full_melt_results_dir`: an existing full-melt run to use as the starting
  state, unless `params.rerun_full_melt=True`.

For detailed parameter descriptions, see [`PARAMS.md`](PARAMS.md). For a visual
overview of the run flow, see [`WORKFLOW.md`](WORKFLOW.md). For module roles and
call order, see [`ARCHITECTURE.md`](ARCHITECTURE.md).

## Replot An Existing Run

Set `just_plots=True` and point `plot_results_dir` to an existing partial-melt
results directory:

```python
run_partial_melt(
    params=PartialMeltParams(),
    run_name="unused_when_just_plotting",
    just_plots=True,
    plot_results_dir="results_partial/apr02/test4_partial_melt",
)
```

## Main Files

- [`run_partial_melt.py`](run_partial_melt.py): public entry point for running
  or replotting a partial-melt chain.
- [`partial_melt_orchestrator.py`](partial_melt_orchestrator.py): orchestration
  layer and home of `PartialMeltParams`.
- [`partial_melt_science.py`](partial_melt_science.py): melt-fraction schedule,
  frozen-core bookkeeping, and melt/solid state transforms.
- [`partial_melt_plot_and_filter_results.py`](partial_melt_plot_and_filter_results.py):
  result table rebuilding, filtering, and plotting.
