# Running GCE Parameter Studies

This directory contains the standard full-melt GCE entry point. Use it for
regular atmosphere-magma ocean-metal equilibrium runs.

For the staged partial-melt workflow, use
[`run_partial_melt/README.md`](run_partial_melt/README.md).

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
PYTHONPATH=. python3 Example/run_gce.py
```

Before running, edit [`run_gce.py`](run_gce.py). Most runs only need:

- `params`: the `GCEParams(...)` science setup and sweep values.
- `run_name`: the label used for the results directory.
- `version`: the solver directory, such as `Sulfur_Version`, `Carbon_Version`,
  or `Sulfur_Nitrogen_Version`.

For detailed parameter descriptions, see [`PARAMS.md`](PARAMS.md). For a visual
overview of the pipeline, see [`WORKFLOW.md`](WORKFLOW.md).

## Replot An Existing Run

Set `just_plots=True` and point `plot_results_dir` to an existing results
directory:

```python
run_gce(
    run_name="my_run",
    params=GCEParams(...),
    version="Sulfur_Version",
    just_plots=True,
    plot_results_dir="results/may04/my_run_sulfur",
)
```

## Version Selection

`version` must match an existing solver directory. The pipeline checks that
directory with `make` before non-plotting runs, then copies the selected solver,
`param.dat`, `chem_input.dat`, and `Gibbs.dat` into each generated case folder.

Sulfur is enabled when the version name contains `Sulfur`; nitrogen is enabled
when it contains `Nitrogen`.

## Main Files

- [`run_gce.py`](run_gce.py): public entry point for running or replotting a
  full GCE study.
- [`gce_orchestrator.py`](gce_orchestrator.py): orchestration layer for build,
  input generation, solving, result collection, and plotting.
- [`PARAMS.md`](PARAMS.md): detailed guide to `run_gce(...)` controls and
  `GCEParams(...)` fields.
- [`WORKFLOW.md`](WORKFLOW.md): compact workflow diagram.
