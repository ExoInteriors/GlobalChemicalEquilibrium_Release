"""Phase and element fraction plotting functions."""
import numpy as np

from .helpers.plot_constants import PHASE_COLORS
from .helpers.plotting_helpers import axis_label, axis_series, draw_stacked_fractions, plot_panels_or_single


# ---------------------------------------------------------------------------
# Phase Mole Fraction Plotting
# ---------------------------------------------------------------------------

def _draw_phase_subplot(ax, subset, axis_key):
    """Draw a stacked phase mole fraction plot on a single axes."""
    phase_cols = ['Moles_atm', 'Moles_silicate', 'Moles_metal']
    if not set(phase_cols) <= set(subset.columns):
        return False
    phase_data = subset[phase_cols]
    valid_mask = phase_data.notna().all(axis=1).to_numpy()
    if not np.any(valid_mask):
        return False
    phase_df = subset.loc[valid_mask].reset_index(drop=True)
    x_vals = axis_series(phase_df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(phase_df))
    M_atm = phase_df['Moles_atm']
    M_sil = phase_df['Moles_silicate']
    M_met = phase_df['Moles_metal']
    total_moles = M_atm + M_sil + M_met
    total_moles = total_moles.replace(0, float('nan'))
    frac_atm = M_atm / total_moles
    frac_sil = M_sil / total_moles
    frac_met = M_met / total_moles
    ax.stackplot(
        x_vals,
        frac_met,
        frac_sil,
        frac_atm,
        labels=['metal', 'silicate', 'atm'],
        colors=[PHASE_COLORS['metal'], PHASE_COLORS['silicate'], PHASE_COLORS['atm']],
    )
    ax.set_ylim(0, 1.0)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel('Phase mole fraction')
    ax.set_title('Phase mole fraction distribution')
    ax.legend(loc='upper left')
    return True


def plot_phase_mole_fractions(df, path, axis_key):
    """Stack metal/silicate/atm phase mole fractions along the chosen axis."""
    if df is None or len(df) == 0:
        return

    def _title(axis_key, value):
        if value is None:
            return "Phase mole fraction distribution"
        # Format panel title based on axis type
        return ("no accreted water" if value == 0 else f"accreted water = {value:.3g}") if axis_key == "HHe" \
            else (f"HHe = {value:.3g}" if axis_key == "Water" else str(value))

    plot_panels_or_single(df, axis_key, _draw_phase_subplot, _title, path, "phase_mole_fractions")


# ---------------------------------------------------------------------------
# Element Fraction Plotting
# ---------------------------------------------------------------------------

def _sulfur_annotations(ax, x_vals, fractions, labels):
    """Add sulfur fraction text labels above the plot."""
    if 'nS' not in labels:
        return
    s_index = labels.index('nS')
    sulfur_fraction = fractions[:, s_index]
    for xi, frac in zip(x_vals, sulfur_fraction):
        if np.isnan(frac):
            continue
        ax.text(xi, 1.01, f"{frac:.3f}", ha='center', va='bottom', fontsize='x-small', color='black')


def _plot_element_panel(ax, subset, axis_key, element_cols):
    """Draw a stacked element fraction plot on a single axes."""
    element_data = subset[element_cols]
    valid_mask = element_data.notna().all(axis=1).to_numpy()
    if not np.any(valid_mask):
        return False
    elem_df = subset.loc[valid_mask].reset_index(drop=True)
    x_vals = axis_series(elem_df, axis_key)
    if len(x_vals) == 0:
        return False

    elem_matrix = elem_df[element_cols].to_numpy(dtype=float)
    # Normalize per-case totals so fractions sum to 1
    totals = elem_matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    elem_frac = elem_matrix / totals

    draw_stacked_fractions(
        ax, x_vals, elem_frac, element_cols,
        use_alpha_cycle=True,
        annotation_fn=_sulfur_annotations,
    )
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel('Element fraction (moles)')
    ax.set_title('Element distribution (normalized per case)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    return True


def plot_element_fractions(df, path, axis_key):
    """Plot mole fractions (normalized to out of 1) of the tracked elements along an axis."""
    element_cols = [c for c in ['nSi', 'nMg', 'nFe', 'nO', 'nH', 'nNa', 'nC', 'nS'] if c in df.columns]
    if len(element_cols) < 3:
        return

    # Closure to pass element_cols to the draw function
    def _draw_with_cols(ax, subset, axis_key):
        return _plot_element_panel(ax, subset, axis_key, element_cols)

    def _title(axis_key, value):
        if value is None:
            return "Element distribution (normalized per case)"
        # Format panel title based on axis type
        return ("no accreted water" if value == 0 else f"accreted water = {value:.3g}") if axis_key == "HHe" \
            else (f"HHe = {value:.3g}" if axis_key == "Water" else str(value))

    plot_panels_or_single(
        df,
        axis_key,
        _draw_with_cols,
        _title,
        path,
        "element_fractions",
    )
