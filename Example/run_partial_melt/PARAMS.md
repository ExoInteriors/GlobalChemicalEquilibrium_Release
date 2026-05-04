# Partial-Melt Parameter Guide

This guide describes the controls used by
[`run_partial_melt.py`](run_partial_melt.py). The user-facing science settings
live in `PartialMeltParams`, which is defined in
[`partial_melt_organizer.py`](partial_melt_organizer.py).

## Typical Edit Pattern

```python
run_partial_melt(
    params=PartialMeltParams(
        f_melt_stop=0.05,
        f_melt_step=0.1,
        rerun_full_melt=True,
        refractory_gas_to_mantle=True,
    ),
    run_name="test",
    just_plots=False,
)
```

## Run Controls

These parameters belong to `run_partial_melt(...)`, not to
`PartialMeltParams(...)`.

### `run_name`

Name for the partial-melt chain. New results are written under:

```text
results_partial/<date>/<run_name>_partial_melt/
```

The `<date>` folder is numeric `YYYYMMDD`. If `run_name` includes a `/`, the
path before the final `/` is used as the results location and the final path
part is used as the run label:

```python
run_name="my_partial_results/test"       # my_partial_results/YYYYMMDD/test_partial_melt/
run_name="/path/to/results_partial/test" # /path/to/results_partial/YYYYMMDD/test_partial_melt/
```

### `just_plots`

Set `False` for a normal run. Set `True` to skip solving and rebuild plots from
an existing partial-melt results directory.

### `verbose`

Set `False` to hide low-level solver and `Gibbs.py` output during the upstream
full-melt run and each partial-melt step. The default is `True`, which preserves
the full solver output.

### `plot_results_dir`

Existing partial-melt results directory used when `just_plots=True`. Example:

```python
plot_results_dir="results_partial/20260504/test4_partial_melt"
```

### `full_melt_results_dir`

Existing full-melt GCE results directory to use as the starting state when
`params.rerun_full_melt=False`.

### `version_full_melt`

Full-melt solver version used when `params.rerun_full_melt=True`. The current
default is `Sulfur_Nitrogen_Version`.

### `axis_list`

Optional explicit plot-axis list. Leave as `None` for the usual partial-melt
axis behavior, which defaults to `["f_melt"]`.

## `PartialMeltParams`

### `full_melt_params`

`GCEParams(...)` object used to create the upstream full-melt case when
`rerun_full_melt=True`. See [`../PARAMS.md`](../PARAMS.md) for the full GCE
parameter guide.

### `f_melt_stop`

Lowest target melt fraction in the chain. Values are fractions, so `0.05` means
stop at 5 percent melt.

Default: `0.05`.

### `f_melt_step`

Step size between target melt fractions. For example, `0.1` creates targets such
as `1.0`, `0.9`, `0.8`, and so on until `f_melt_stop`.

Use a smaller value if you need specific near-initial or near-final points, such
as `0.99` or `0.01`.

Default: `0.05`.

### `partial_start_t_sme`

Temperature used for the partial-melt silicate-metal equilibrium state.

Default: `2000.0`.

### `refractory_gas_to_mantle`

When `False`, the standard partial-melt gas network remains active.

When `True`, the workflow builds the partial-melt solver with
`No_Refractory_Gas_Version/Equations.py`, excluding refractory gas species such
as `Fe_gas`, `Mg_gas`, `Na_gas`, `SiO_gas`, and `SiH4_gas` from the active
atmospheric network.

Default: `False`.

### `build`

When `True`, rebuild the partial-melt solver before running the chain. This is
especially important when switching `refractory_gas_to_mantle`.

Default: `True`.

### `volatile_retention_in_solid`

Controls whether volatile silicate species freeze into the solid reservoir.

When `False`, volatile silicate species stay in the melt during ordinary
partial-melt steps. When `True`, they are allowed to enter the solid reservoir
with the freezing material.

Default: `False`.

### `freeze_solid`

Controls whether the active melt reservoir shrinks through the chain.

When `True`, each step freezes additional silicate according to the target melt
fraction. When `False`, the workflow still solves chained states but does not
remove additional active silicate mass.

Default: `True`.

### `rerun_full_melt`

When `True`, run a fresh upstream full-melt GCE case using `full_melt_params`
and `version_full_melt`.

When `False`, load an existing full-melt run from `full_melt_results_dir`.

Default: `False`.

## Freeze Bookkeeping

The workflow uses the original post-core-freeze silicate melt mass as the
reference mass, stored as `M_silicate_ref`. That reference is not redefined after
melt-gas rebalancing.

For example, if `M_silicate_ref = 100` and the target moves from `f_melt=1.0` to
`f_melt=0.95`, the workflow freezes 5 units of the original reference mass. If
rebalancing changes the active silicate mass to 94 or 89 afterward, the next
target still uses the same original reference mass to determine cumulative
freezing.

This is why a requested target can sometimes be skipped: the current bookkeeping
state may already be at or beyond that melt fraction.
