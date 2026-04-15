"""Functions that help with science processing when creating plots"""

import numpy as np
import pandas as pd

from tools.calc_fO2 import get_delta_IW, log10_fO2_IW_hirschmann2021

from .plot_constants import PHASE_MOLES_COLUMNS, PHASE_ORDER, SULFUR_SPECIES_MW, format_species_label

def get_delta_iw_series(df):
    """Return ΔIW (log fO2 model − log fO2_IW) series for plotting."""
    if df is None or df.empty:
        return np.array([])
    required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
    if not required.issubset(df.columns):
        return np.full(len(df), np.nan, dtype=float)

    temperature = df["T_SME"].to_numpy(dtype=float)
    n_melt = df["Moles_silicate"].to_numpy(dtype=float)
    n_metal = df["Moles_metal"].to_numpy(dtype=float)
    n_feo = df["FeO_silicate"].to_numpy(dtype=float) * n_melt
    n_fe = df["Fe_metal"].to_numpy(dtype=float) * n_metal

    # Pressure choice is arbitrary here because the IW buffer terms cancel in ΔIW.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        delta_iw = get_delta_IW(10.0, temperature, n_feo, n_melt, n_fe, n_metal)
    return np.where(np.isfinite(delta_iw), delta_iw, np.nan)


def get_log10_fO2_series(df):
    """Return log10(fO2) in bar for each row."""
    if df is None or df.empty:
        return np.array([])
    required = {"T_SME", "P_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
    if not required.issubset(df.columns):
        return np.full(len(df), np.nan, dtype=float)

    temperature = df["T_SME"].to_numpy(dtype=float)
    pressure_gpa = df["P_SME"].to_numpy(dtype=float)
    n_melt = df["Moles_silicate"].to_numpy(dtype=float)
    n_metal = df["Moles_metal"].to_numpy(dtype=float)
    n_feo = df["FeO_silicate"].to_numpy(dtype=float) * n_melt
    n_fe = df["Fe_metal"].to_numpy(dtype=float) * n_metal

    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        a_feo = n_feo / n_melt
        a_fe = n_fe / n_metal
        log10_fO2_iw = log10_fO2_IW_hirschmann2021(pressure_gpa, temperature)
        log10_fO2 = log10_fO2_iw + 2.0 * np.log10(a_feo / a_fe)
    return np.where(np.isfinite(log10_fO2), log10_fO2, np.nan)


def get_matm_mplanet_series(df):
    """Return the atmospheric mass fraction: Matm / Mtotal."""
    if df is None or df.empty:
        return np.array([])
    from .plotting_helpers import load_atomic_weights, mass_arrays

    atomic_weights = load_atomic_weights()
    _, _, _, grams_atm, _, _, total_mass = mass_arrays(df, atomic_weights)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(total_mass > 0, grams_atm / total_mass, 0.0)
    return np.where(np.isfinite(ratio), ratio, 0.0)


def get_f_solid_series(df):
    """Return the actual silicate solid fraction in percent for plotting."""
    if df is None or df.empty:
        return np.array([])
    if "M_frozen_solid" in df.columns and "Moles_silicate" in df.columns:
        from .plotting_helpers import load_atomic_weights, mass_arrays

        atomic_weights = load_atomic_weights()
        _, _, _, _, grams_silicate, _, _ = mass_arrays(df, atomic_weights)
        grams_solid = np.nan_to_num(
            df["M_frozen_solid"].to_numpy(dtype=float),
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        )
        if "M_silicate_ref" in df.columns:
            reference_mass = pd.to_numeric(df["M_silicate_ref"], errors="coerce").to_numpy(dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                f_solid = np.where(reference_mass > 0, grams_solid / reference_mass, 0.0)
            return np.where(np.isfinite(f_solid), 100.0 * f_solid, np.nan)

        total_silicate = grams_solid + grams_silicate
        reference_mass = np.nan
        if "f_melt" in df.columns:
            f_melt = pd.to_numeric(df["f_melt"], errors="coerce").to_numpy(dtype=float)
            finite_mask = np.isfinite(f_melt) & np.isfinite(total_silicate)
            if np.any(finite_mask):
                reference_mass = total_silicate[finite_mask][np.argmax(f_melt[finite_mask])]
        if not np.isfinite(reference_mass) or reference_mass <= 0.0:
            finite_total = total_silicate[np.isfinite(total_silicate) & (total_silicate > 0)]
            if finite_total.size > 0:
                reference_mass = np.max(finite_total)

        with np.errstate(divide="ignore", invalid="ignore"):
            f_solid = np.where(reference_mass > 0, grams_solid / reference_mass, 0.0)
        return np.where(np.isfinite(f_solid), 100.0 * f_solid, np.nan)

    if "f_melt" not in df.columns:
        return np.zeros(len(df))
    f_melt = df["f_melt"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        f_solid = 1.0 - f_melt
    return np.where(np.isfinite(f_solid), 100.0 * f_solid, np.nan)




def prepare_atmosphere_partial_pressures(subset, columns, axis_key):
    """Return sorted x values and species partial pressures for atmosphere plots."""
    if subset is None or subset.empty:
        return None
    from .plotting_helpers import axis_series

    gas_df = subset.reindex(columns=columns, fill_value=0.0).astype(float).fillna(0.0)
    p_total = pd.to_numeric(subset.get("Pstd", 0.0), errors="coerce").to_numpy(dtype=float)
    partial_pressures = gas_df.mul(p_total, axis=0)
    if np.nan_to_num(partial_pressures.to_numpy(), nan=0.0).sum() == 0:
        return None

    x_vals = axis_series(subset, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(subset))

    valid_mask = np.isfinite(x_vals) & np.isfinite(p_total)
    if not np.any(valid_mask):
        return None

    x_vals = np.asarray(x_vals)[valid_mask]
    partial_pressures = partial_pressures.loc[valid_mask].reset_index(drop=True)
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    partial_pressures = partial_pressures.iloc[order].reset_index(drop=True)
    labels = [format_species_label(col) for col in columns]
    return x_vals, partial_pressures.to_numpy(dtype=float), labels


def detect_matm_dual_axis(df, axis_key):
    """Detect dual x-axis configuration for Matm/Mplanet plots."""
    if axis_key != "Matm_Mplanet":
        return None
    from .plotting_helpers import axis_series

    hhe_vals = axis_series(df, "HHe")
    water_vals = axis_series(df, "Water")
    hhe_varied = len(np.unique(hhe_vals[np.isfinite(hhe_vals)])) > 1 if hhe_vals.size > 0 else False
    water_varied = len(np.unique(water_vals[np.isfinite(water_vals)])) > 1 if water_vals.size > 0 else False

    if hhe_varied:
        return {
            "bottom_axis_key": "HHe",
            "label": r"Accreted H from primordial gas (wt \%)",
            "bottom_vals": hhe_vals,
            "matm_vals": get_matm_mplanet_series(df),
        }
    if water_varied:
        return {
            "bottom_axis_key": "Water",
            "label": r"Accreted water after formation (wt \%)",
            "bottom_vals": water_vals,
            "matm_vals": None,
        }
    return None


AXIS_CONFIG: dict[str, dict] = {
    "HHe": {"axis_key": "HHe", "label": r"Accreted H from primordial gas (wt \%)"},
    "P_GPa": {"axis_key": "P_SME", "label": "SME pressure (GPa)"},
    "delta_IW": {"column": "delta_IW", "label": r"$\Delta$IW (log $f_{\mathrm{O}_2}$ model $-$ log $f_{\mathrm{O}_2,\mathrm{IW}}$)"},
    "log10_fO2": {"column": "log10_fO2", "label": r"$\log_{10}(f_{\mathrm{O}_2})$ (bar)"},
    "Matm_Mplanet": {"getter": get_matm_mplanet_series, "label": "Atmosphere mass fraction / planet mass"},
}

SLICE_CONFIG: dict[str, dict] = {
    "Planetmass": {
        "axis_key": "Planetmass",
        "label": "Planet mass",
        "unit": r"$M_\oplus$",
        "format": lambda value: rf"M={value:.1f} $M_\oplus$",
    },
    "P_SME": {
        "axis_key": "P_SME",
        "label": "SME pressure",
        "unit": "GPa",
        "format": lambda value: f"P={value:.1f} GPa",
    },
    "HHe": {
        "axis_key": "HHe",
        "label": "Accreted H",
        "unit": r"wt \%",
        "format": lambda value: rf"H={value:.3g} wt\%",
    },
}


def get_config_values(df, config: dict) -> np.ndarray:
    """Extract values from a DataFrame using an AXIS_CONFIG or SLICE_CONFIG entry."""
    from .plotting_helpers import axis_series

    if "getter" in config:
        return config["getter"](df)
    if "column" in config:
        return df[config["column"]].to_numpy(dtype=float)
    if "axis_key" in config:
        return axis_series(df, config["axis_key"])
    raise ValueError("Config must have 'getter', 'column', or 'axis_key'")


def accumulate_element_by_phase(df, element, weights=None, phase_moles=None):
    """Accumulate element contributions by phase from species data."""
    from .plot_constants import ELEMENT_SPECIES

    length = len(df)
    result = {
        "atm": np.zeros(length, dtype=float),
        "silicate": np.zeros(length, dtype=float),
        "metal": np.zeros(length, dtype=float),
    }

    for species, coeff in ELEMENT_SPECIES.get(element, []):
        arr = df[species].to_numpy(dtype=float) if species in df.columns else np.zeros(length, dtype=float)
        phase = "atm" if species.endswith("_gas") else "metal" if species.endswith("_metal") else "silicate"
        weight = weights.get(species, 1.0) if weights else 1.0
        if weight == 0.0:
            continue
        phase_mult = phase_moles.get(phase, 1.0) if phase_moles else 1.0
        result[phase] = result[phase] + coeff * arr * weight * phase_mult

    return result


def sulfur_phase_mass_fractions(df):
    """Return sulfur mass fractions per phase (metal/silicate/atm)."""
    if df is None or df.empty:
        empty = np.array([], dtype=float)
        return {"metal": empty, "silicate": empty, "atm": empty}

    phase_moles = {
        "atm": df["Moles_atm"].to_numpy(dtype=float) if "Moles_atm" in df.columns else np.zeros(len(df), dtype=float),
        "silicate": df["Moles_silicate"].to_numpy(dtype=float) if "Moles_silicate" in df.columns else np.zeros(len(df), dtype=float),
        "metal": df["Moles_metal"].to_numpy(dtype=float) if "Moles_metal" in df.columns else np.zeros(len(df), dtype=float),
    }
    masses = accumulate_element_by_phase(df=df, element="S", weights=SULFUR_SPECIES_MW, phase_moles=phase_moles)
    total = masses["atm"] + masses["silicate"] + masses["metal"]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_atm = np.where(total > 0, masses["atm"] / total, 0.0)
        frac_sil = np.where(total > 0, masses["silicate"] / total, 0.0)
        frac_met = np.where(total > 0, masses["metal"] / total, 0.0)
    return {"metal": frac_met, "silicate": frac_sil, "atm": frac_atm}


def compute_phase_mole_fractions(df, axis_key):
    """Compute phase mole fractions sorted by the given axis key."""
    if df is None or df.empty:
        return None
    from .plotting_helpers import axis_series

    phase_cols = ["Moles_atm", "Moles_silicate", "Moles_metal"]
    if not set(phase_cols) <= set(df.columns):
        return None
    phase_data = df[phase_cols]
    valid_mask = phase_data.notna().all(axis=1).to_numpy()
    if not np.any(valid_mask):
        return None
    phase_df = df.loc[valid_mask].reset_index(drop=True)

    x_vals = axis_series(phase_df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(phase_df))

    m_atm = phase_df["Moles_atm"].to_numpy(dtype=float)
    m_sil = phase_df["Moles_silicate"].to_numpy(dtype=float)
    m_met = phase_df["Moles_metal"].to_numpy(dtype=float)
    total_moles = m_atm + m_sil + m_met

    finite_mask = np.isfinite(total_moles) & (total_moles > 0)
    if not np.any(finite_mask):
        return None
    x_vals = np.asarray(x_vals)[finite_mask]
    m_atm = m_atm[finite_mask]
    m_sil = m_sil[finite_mask]
    m_met = m_met[finite_mask]
    total_moles = total_moles[finite_mask]

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    frac_atm = (m_atm / total_moles)[order]
    frac_sil = (m_sil / total_moles)[order]
    frac_met = (m_met / total_moles)[order]
    return x_vals, frac_atm, frac_sil, frac_met


def compute_phase_mass_fractions(df, axis_key):
    """Compute phase mass fractions sorted by the given axis key."""
    if df is None or df.empty:
        return None
    from .plotting_helpers import axis_series, load_atomic_weights, mass_arrays

    mu = load_atomic_weights()
    _, _, _, grams_atm, grams_silicate, grams_metal, total_mass = mass_arrays(df, mu)
    if not np.any(total_mass > 0):
        return None

    x_vals = axis_series(df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(df))

    valid_mask = np.isfinite(x_vals) & np.isfinite(total_mass) & (total_mass > 0)
    if not np.any(valid_mask):
        return None

    x_vals = np.asarray(x_vals)[valid_mask]
    grams_atm = grams_atm[valid_mask]
    grams_silicate = grams_silicate[valid_mask]
    grams_metal = grams_metal[valid_mask]
    total_mass = total_mass[valid_mask]

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    grams_atm = grams_atm[order]
    grams_silicate = grams_silicate[order]
    grams_metal = grams_metal[order]
    total_mass = total_mass[order]

    total_safe = np.where(total_mass == 0, np.nan, total_mass)
    frac_atm = np.nan_to_num(grams_atm / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
    frac_silicate = np.nan_to_num(grams_silicate / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
    frac_metal = np.nan_to_num(grams_metal / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
    return x_vals, frac_atm, frac_silicate, frac_metal


def _atomic_weight_for_element_column(mu: dict, col: str) -> float:
    """Resolve atomic weight (g/mol) for an inventory column like ``nSi``."""
    element = col[1:]
    aw = mu.get(f"{element}_metal") or mu.get(f"{element}_gas")
    if aw is None:
        diatomic_mw = mu.get(f"{element}2_gas")
        if diatomic_mw is not None:
            aw = diatomic_mw / 2.0
    if aw is None:
        raise ValueError(f"No atomic weight found for element '{element}' in Molecular_Weight.dat")
    return float(aw)


def element_gce_mass_scores(gce_row, element_cols):
    """Per-element mass at the saved GCE state (moles × atomic weight) for ranking partial-melt plots."""
    if gce_row is None:
        return {}
    from .plotting_helpers import load_atomic_weights

    mu = load_atomic_weights()
    out = {}
    for col in element_cols:
        aw = _atomic_weight_for_element_column(mu, col)
        n = float(pd.to_numeric(gce_row.get(col, np.nan), errors="coerce"))
        out[col] = n * aw if np.isfinite(n) and n > 0.0 else 0.0
    return out


def silicate_species_gce_masses(gce_row, columns):
    """Absolute silicate species mass in the melt at the saved GCE state; ``*_solid_frac`` columns are 0."""
    if gce_row is None or not columns:
        return {}
    from .plotting_helpers import load_atomic_weights

    mu = load_atomic_weights()
    moles_melt = float(
        np.nan_to_num(gce_row.get("Moles_silicate", 0.0), nan=0.0, posinf=0.0, neginf=0.0)
    )
    out = {}
    for col in columns:
        if col.endswith("_solid_frac"):
            out[col] = 0.0
            continue
        frac = float(np.nan_to_num(gce_row.get(col, np.nan), nan=0.0))
        out[col] = moles_melt * frac * mu.get(col, 0.0)
    return out


def compute_element_weight_fractions(df, element_cols, axis_key, *, sort_by_mean: bool = True):
    """Compute element weight fractions along ``axis_key``, sorted by axis and optionally by mean wt%."""
    if df is None or df.empty:
        return None
    from .plotting_helpers import axis_series, load_atomic_weights

    mu = load_atomic_weights()
    element_data = df[element_cols]
    valid_mask = element_data.notna().all(axis=1).to_numpy()
    if not np.any(valid_mask):
        return None
    elem_df = df.loc[valid_mask].reset_index(drop=True)

    x_vals = axis_series(elem_df, axis_key)
    if len(x_vals) == 0:
        return None

    elem_matrix = elem_df[element_cols].to_numpy(dtype=float)
    mass_matrix = np.zeros_like(elem_matrix)
    for i, col in enumerate(element_cols):
        aw = _atomic_weight_for_element_column(mu, col)
        mass_matrix[:, i] = elem_matrix[:, i] * aw

    totals = mass_matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    wt_frac = mass_matrix / totals

    if sort_by_mean:
        avg_wt = wt_frac.mean(axis=0)
        sort_idx = np.argsort(-avg_wt)
    else:
        sort_idx = np.arange(len(element_cols), dtype=int)
    sorted_frac = wt_frac[:, sort_idx]
    sorted_labels = [element_cols[i][1:] for i in sort_idx]

    x_vals = np.asarray(x_vals)
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    sorted_frac = sorted_frac[order]
    return x_vals, sorted_frac, sorted_labels


def prepare_phase_fractions(df, columns, axis_key):
    """Return sorted x values, mass fractions, and labels for plotting."""
    if df is None or df.empty:
        return None
    from .plotting_helpers import axis_series, load_atomic_weights

    phase_df = df.reindex(columns=columns, fill_value=0.0).astype(float).fillna(0.0)
    mu = load_atomic_weights()

    def _weight_key(name: str) -> str:
        return name[:-len("_solid_frac")] if name.endswith("_solid_frac") else name

    weights = {_weight_key(name): mu.get(_weight_key(name), 0.0) for name in columns if name in phase_df.columns}
    missing_weights = [col for col in columns if _weight_key(col) not in mu]
    if missing_weights:
        raise ValueError(f"Missing molar masses for: {', '.join(missing_weights)}")
    weighted_df = phase_df.mul([weights.get(_weight_key(col), 0.0) for col in phase_df.columns], axis=1)
    total_mass = weighted_df.sum(axis=1).to_numpy()
    total_mass = np.nan_to_num(total_mass, nan=0.0, posinf=0.0, neginf=0.0)
    if not np.any(total_mass > 0):
        return None

    x_vals = axis_series(df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(df))

    valid_mask = np.isfinite(x_vals) & np.isfinite(total_mass)
    if not np.any(valid_mask):
        return None
    x_vals = np.asarray(x_vals)[valid_mask]
    weighted_df = weighted_df.loc[valid_mask].reset_index(drop=True)
    total_mass = total_mass[valid_mask]

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    weighted_df = weighted_df.iloc[order].reset_index(drop=True)
    total_mass = total_mass[order]

    total_safe = np.where(total_mass == 0, np.nan, total_mass)
    fractions = weighted_df.to_numpy() / total_safe[:, None]
    fractions = np.nan_to_num(fractions, nan=0.0, posinf=0.0, neginf=0.0)
    labels = [format_species_label(col) for col in columns]
    return x_vals, fractions, labels


def prepare_mole_fractions(df, columns, axis_key):
    """Return sorted x values, raw mole fractions from results.dat, and labels."""
    if df is None or df.empty:
        return None
    from .plotting_helpers import axis_series

    phase_df = df.reindex(columns=columns, fill_value=0.0).astype(float).fillna(0.0)
    if phase_df.sum().sum() == 0:
        return None

    x_vals = axis_series(df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(df))

    valid_mask = np.isfinite(x_vals)
    if not np.any(valid_mask):
        return None

    x_vals = np.asarray(x_vals)[valid_mask]
    phase_df = phase_df.loc[valid_mask].reset_index(drop=True)

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    phase_df = phase_df.iloc[order].reset_index(drop=True)
    labels = [format_species_label(col) for col in columns]
    return x_vals, phase_df.to_numpy(), labels


def compute_and_filter(df, series_fn, column_name, required_cols, label):
    """Compute a derived series, filter to valid rows, and return sorted DataFrame."""
    values = series_fn(df)
    if np.all(np.isnan(values)):
        missing = sorted(required_cols - set(df.columns))
        if missing:
            raise ValueError(f"Cannot compute {label} because the following columns are missing: {', '.join(missing)}")

    n_melt = df["Moles_silicate"].to_numpy(dtype=float) if "Moles_silicate" in df.columns else np.zeros(len(df))
    valid = np.isfinite(values) & np.isfinite(n_melt) & (n_melt > 0)
    if not np.any(valid):
        raise ValueError(f"No valid {label} values could be computed for the dataset.")

    subset = df.loc[valid].copy()
    subset[column_name] = values[valid]
    return subset.sort_values(column_name).reset_index(drop=True)


def compute_atm_co_ratio(df):
    """Return the atmospheric C/O mole ratio."""
    if "upperCO" in df.columns:
        return df["upperCO"].to_numpy(dtype=float)

    n = len(df)
    c_atm = np.zeros(n, dtype=float)
    for col in ("CH4_gas", "CO2_gas", "CO_gas"):
        if col in df.columns:
            c_atm += df[col].to_numpy(dtype=float)

    o_atm = np.zeros(n, dtype=float)
    if "CO_gas" in df.columns:
        o_atm += df["CO_gas"].to_numpy(dtype=float)
    if "CO2_gas" in df.columns:
        o_atm += 2.0 * df["CO2_gas"].to_numpy(dtype=float)
    if "O2_gas" in df.columns:
        o_atm += 2.0 * df["O2_gas"].to_numpy(dtype=float)
    if "H2O_gas" in df.columns:
        o_atm += df["H2O_gas"].to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        co_ratio = np.where(o_atm > 0, c_atm / o_atm, np.nan)
    return co_ratio


def sulfur_phase_mole_fractions(df):
    """Return sulfur mole fractions in each phase as a dict {phase: array}."""
    phase_moles = {
        phase: df[col].to_numpy(dtype=float) if col in df.columns else np.zeros(len(df))
        for phase, col in PHASE_MOLES_COLUMNS.items()
    }
    moles = accumulate_element_by_phase(df, "S", phase_moles=phase_moles)
    total = sum(moles[phase] for phase in PHASE_ORDER)
    total = np.where(total == 0, np.nan, total)
    return {phase: moles[phase] / total for phase in PHASE_ORDER}
