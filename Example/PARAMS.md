# GCE Parameter Guide

This guide describes the user-facing parameters for the standard full-melt GCE
entry point:

```bash
PYTHONPATH=. python3 Example/run_gce.py
```

Most science runs are controlled by the `GCEParams(...)` object passed to
`run_gce()` in [`run_gce.py`](run_gce.py). Each array-valued field can contain
one value or several values. When several fields contain multiple values, the
pipeline runs the Cartesian product of those choices, so the number of cases can
grow quickly.

For example, 6 water values and 6 S/H values produce 36 solver cases before any
other multi-valued parameters are added.

## Typical Edit Pattern

```python
run_gce(
    run_name="my_run",
    params=GCEParams(
        Planetmassarray=np.array([5.0]),
        T_AMOI_array=np.array([2500.0]),
        T_SME_array=np.array([3000.0]),
        tarWaterarray=np.array([0.0, 0.02, 0.04]),
        tarHHearray=np.array([0.03]),
        tarDiskCOarray=np.array([0.5]),
        tarDiskSHarray=np.array([1.335e-5]),
    ),
    version="Sulfur_Version",
    just_plots=False,
)
```

Use one-element arrays for fixed values and longer arrays for parameter sweeps.

## Run Controls

These parameters belong to `run_gce(...)`, not to `GCEParams(...)`.

### `run_name`

Human-readable name for the run. New results are written under:

```text
results/<date>/<run_name>_<version-short>/
```

For example, `run_name="testy"` with `version="Sulfur_Version"` produces a
folder like `results/20260504/testy_sulfur/`. If the folder already exists
during a fresh run, the pipeline creates an incremented name such as
`testy_2_sulfur`.

The `<date>` folder is numeric `YYYYMMDD`. If `run_name` includes a `/`, the
path before the final `/` is used as the results location and the final path
part is used as the run label:

```python
run_name="my_results/testy"       # my_results/YYYYMMDD/testy_sulfur/
run_name="/path/to/results/testy" # /path/to/results/YYYYMMDD/testy_sulfur/
```

### `version`

Selects which chemistry/solver directory to build and run. Common choices are:

- `Standard_Version`
- `Carbon_Version`
- `Sulfur_Version`
- `Sulfur_Nitrogen_Version`

The selected version controls which elements and species are available. If the
version name does not include `Sulfur`, sulfur inputs are set to zero. If the
version name does not include `Nitrogen`, nitrogen inputs are set to zero.

### `just_plots`

Set `False` for a normal run: create inputs, copy solver files, solve all cases,
collect results, and generate plots.

Set `True` to skip solver work and regenerate plots from an existing results
directory. When `just_plots=True`, also set `plot_results_dir`.

### `plot_results_dir`

Path to an existing results directory, relative to the repository root, used
when replotting. Example:

```python
plot_results_dir="results/20260504/testy_sulfur"
```

### `only_sulfur_plots`

Controls plot selection. With `True`, the plotting step focuses on the
sulfur-specific plot set. With `False`, it also runs the broader plot suite.

### `axis_list`

Optional explicit list of plot axes. Usually leave this as `None`; the pipeline
infers axes from the `GCEParams` fields that have multiple values. Use this only
when you need to override the automatic choice.

Known axis keys include `HHe`, `Water`, `P_AMOI`, `P_SME`, `T_AMOI`, `T_SME`,
`Planetmass`, `O`, `f_melt`, `delta_IW`, and `Matm_Mplanet`.

## Science Parameters

These fields belong inside `GCEParams(...)`.

### `Planetmassarray`

Planet mass in Earth masses. This value is written to `chem_input.dat` as
`Mplanet_Mearth` and to `planetary_params.dat` as `Mplanet_Mearth`.

Use this when changing the planet size assumed by the run. The input generator
also computes an internal `P_Carbon` estimate from planet mass for masses above
`1.64 Mearth`; that value is written to `planetary_params.dat`.

Default in `GCEParams`: `np.array([5.0])`.

### `FakeMolesTotal`

Normalization scale for the starting bulk composition. It sets the total
synthetic inventory used to convert chondritic abundance presets into moles.

This is not usually a science sweep parameter. It mainly controls the numerical
scale of the generated mole counts before water, H/He gas, and disk volatile
additions are applied.

Default: `10e3`.

### `T_AMOI_array`

Atmosphere-magma ocean interface temperature in Kelvin. It is written to
`chem_input.dat` as `T_AMOI` and to `planetary_params.dat` as `T_surf`.

Use this to sweep or set the surface/interface temperature of the full-melt
equilibrium state.

Default in `GCEParams`: `np.array([2500.0])`.

### `T_SME_array`

Silicate-metal equilibrium temperature in Kelvin. It is written to
`chem_input.dat` as `T_SME` and to `planetary_params.dat` as `T_CMB`.

Use this to control the deeper metal-silicate equilibration temperature.

Default in `GCEParams`: `np.array([3000.0])`.

### `tarmgsiarray`

Target Mg/Si molar ratio. The input generator first computes baseline silicon
moles, then sets:

```text
nMg = tarmgsiarray value * nSi
```

Use this to enrich or deplete magnesium relative to silicon while keeping the
same baseline silicon inventory.

Default in `GCEParams`: `np.array([1.0])`.

### `tarfesiarray`

Target Fe/Si molar ratio. The input generator first computes baseline silicon
moles, then sets:

```text
nFe = tarfesiarray value * nSi
```

Use this to change the metal budget relative to silicon.

Default in `GCEParams`: `np.array([1.0])`.

### `tarWaterarray`

Added water mass fraction. This value represents water added after the baseline
chondritic composition is assembled.

For each value, the code computes a water mass from the initial rocky mass and
adds the corresponding hydrogen and oxygen:

```text
Mwater = f_water * Mrock / (1 - f_water - f_HHe)
nO += moles of H2O
nH += 2 * moles of H2O
```

The code requires:

```text
tarWaterarray value + tarHHearray value < 1
```

because the remaining fraction is the rocky inventory. Values are fractions, not
percentages, so `0.02` means 2 wt%.

Default in `GCEParams`: `np.array([0.0])`.

### `tarHHearray`

Accreted primordial H/He mass fraction. The current implementation uses this as
the mass fraction of primordial gas added after water is included.

The code adds hydrogen from the primordial gas and assumes 2 percent of that gas
mass is carried by C and O. That C/O addition is controlled by `tarDiskCOarray`.

Values are fractions, not percentages, so `0.03` means 3 wt%.

Default in `GCEParams`: `np.array([0.03])`.

### `tarDiskCOarray`

Disk gas C/O ratio used when adding the C and O carried by the primordial gas.
It affects how the fixed 2 percent heavy-element component of the H/He gas is
split between carbon and oxygen.

Higher values add relatively more carbon for the same primordial gas mass; lower
values add relatively more oxygen.

Default in `GCEParams`: `np.array([0.5])`.

### `tarDiskSHarray`

Disk gas S/H ratio. This only matters for versions whose name includes
`Sulfur`. For non-sulfur versions, the input generator uses `0.0` regardless of
this field.

When sulfur is enabled, the code adds sulfur in proportion to the primordial
hydrogen inventory:

```text
nS += tarDiskSHarray value * nH_prim
```

Default in `GCEParams`: `np.array([1.335e-5])`.

### `P_SME_array`

Silicate-metal equilibrium pressure parameter. For non-nitrogen sulfur versions,
this is written to `chem_input.dat` and `parametersAll.dat` as `P_SME`.

In the plotting helpers, `P_SME` is labelled in GPa. Keep units consistent with
the selected solver version and its `param.dat`/equation assumptions.

Default in `GCEParams`: `np.array([10.0])`.

### `UseCondriticComp`

Selects how the baseline chondritic composition is interpreted.

Allowed values:

- `"molar fraction"`
- `"mass fraction"`

With `"molar fraction"`, the code uses the Ed Young/Kallemeyn/Wasson-style
baseline and multiplies each element fraction by `FakeMolesTotal`.

With `"mass fraction"`, the code uses a mass-fraction preset and converts each
element into moles using molecular/atomic weights from `Molecular_Weight.dat`.

Default: `"molar fraction"`.

### `UseCondriticPreset`

Selects the baseline composition preset.

For `UseCondriticComp="molar fraction"`, the active preset is the Ed Young
molar-fraction composition.

For `UseCondriticComp="mass fraction"`, allowed values are:

- `"allegre"`
- `"javoy"`

Default: `"ed_young"`.

If you switch `UseCondriticComp` to `"mass fraction"`, also switch
`UseCondriticPreset` to either `"allegre"` or `"javoy"`; `"ed_young"` is only
valid for the molar-fraction path.

## Version-Dependent Elements

Sulfur and nitrogen are version-gated.

If `version` contains `Sulfur`, sulfur remains active and the code writes `nS`
to `chem_input.dat`. If not, sulfur is set to zero and sulfur-specific sweeps do
not create extra cases.

If `version` contains `Nitrogen`, nitrogen remains active and the code writes
`nN` to `chem_input.dat`. If not, nitrogen is set to zero. The current
`Sulfur_Nitrogen_Version` path does not write `P_SME` as a solver parameter.

## Currently Inactive Or Legacy Fields

These fields exist on `GCEParams`, but the standard `create.py` path currently
does not use them to change generated `chem_input.dat` files:

- `P_AMOI_array`
- `tarFearray`
- `tarOarray`

They may appear in plotting metadata or older workflows, but changing them alone
will not change a standard full-melt solver input in the current pipeline.

## Output Files To Check

After input generation, each case folder contains:

- `chem_input.dat`: the actual solver-facing chemistry and thermodynamic input.
- `planetary_params.dat`: a human-readable record of the generated planetary
  parameters.

The run directory also contains:

- `parametersAll.dat`: compact numerical table consumed by downstream steps.
- `summary_chem_input_GEC.csv`: one row per generated case, including the
  parameter choices and any input-generation failure messages.

When in doubt, inspect `summary_chem_input_GEC.csv` first. It is the fastest way
to confirm that your parameter sweep expanded the way you expected.
