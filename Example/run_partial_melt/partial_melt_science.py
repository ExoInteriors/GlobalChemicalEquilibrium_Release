from __future__ import annotations

"""Science helpers for the partial-melt, including calculating solid, amount frozen, and moving specific species around.
"""

import numpy as np
import pandas as pd

from Example.plots.helpers.plotting_helpers import load_atomic_weights
from Example.plots.helpers.plot_constants import ELEMENT_SPECIES, ELEMENTS, METAL_COLUMNS, SILICATE_COLUMNS, \
    VOLATILE_SILICATE_COLUMNS


# Melt-fraction schedule helpers
# These decide which target f_melt values the chained run will attempt.
def generate_f_melt_schedule(
    f_melt_stop: float,
    f_melt_step: float,
) -> list[float]:
    """Return the inclusive melt-fraction schedule from 1.0 down to the stop value."""
    if not 0.0 <= f_melt_stop <= 1.0:
        raise ValueError(f"f_melt_stop must be between 0 and 1, got {f_melt_stop}")
    if not 0.0 < f_melt_step <= 1.0:
        raise ValueError(f"f_melt_step must be > 0 and <= 1, got {f_melt_step}")
    schedule = [round(value, 12) for value in np.arange(1.0, f_melt_stop, -f_melt_step)]
    if not schedule or round(f_melt_stop, 12) != round(schedule[-1], 12):
        schedule.append(round(float(f_melt_stop), 12))
    # Deduplicate and sort descending to keep the schedule stable across float roundoff.
    return sorted({round(value, 12) for value in schedule}, reverse=True)


# State-mass summaries
# These read one solved state and summarize the active silicate reservoir.
def compute_active_silicate_mass(row: pd.Series) -> float:
    """Compute active silicate mass from one solved state row."""
    mu = load_atomic_weights()
    moles_silicate = float(pd.to_numeric(row.get("Moles_silicate", 0.0), errors="coerce"))
    grams_per_mole_silicate = 0.0
    for species in SILICATE_COLUMNS:
        if species not in row.index:
            continue
        species_fraction = float(pd.to_numeric(row.get(species, 0.0), errors="coerce"))
        grams_per_mole_silicate += species_fraction * mu.get(species, 0.0)
    return moles_silicate * grams_per_mole_silicate


def compute_remaining_melt_fraction_from_bookkeeping(
    row: pd.Series,
    reference_mass: float,
) -> float:
    """Return remaining melt fraction from frozen-mass bookkeeping.

    This is the physically consistent schedule quantity for the chained
    partial-melt workflow: how much of the original silicate reference mass has
    not yet been frozen into the cumulative solid reservoir.

    We intentionally do not infer this from the post-equilibrium silicate phase
    mass, because volatile/interior mass can repartition into the atmosphere
    during each rebalance step without representing additional freezing.
    """
    if reference_mass <= 0.0:
        return 0.0

    frozen_mass = float(pd.to_numeric(row.get("M_frozen_solid", 0.0), errors="coerce"))
    if np.isfinite(frozen_mass):
        remaining = 1.0 - (frozen_mass / reference_mass)
        return float(np.clip(remaining, 0.0, 1.0))

    active_mass = compute_active_silicate_mass(row)
    return float(np.clip(active_mass / reference_mass, 0.0, 1.0))


def compute_additional_freeze_mass_from_bookkeeping(
    row: pd.Series,
    reference_mass: float,
    target_f_melt: float,
) -> float:
    """Return the additional solid mass that should freeze in this step.

    The schedule is defined relative to the original reference melt mass, so a
    target ``f_melt`` means the cumulative frozen mass should become

    ``(1 - target_f_melt) * reference_mass``.

    The incremental freeze for this step is therefore the target cumulative
    frozen mass minus the already-frozen cumulative mass.
    """
    if reference_mass <= 0.0:
        return 0.0

    target_cumulative_frozen = (1.0 - float(target_f_melt)) * reference_mass
    current_cumulative_frozen = float(pd.to_numeric(row.get("M_frozen_solid", 0.0), errors="coerce"))
    additional_freeze = target_cumulative_frozen - current_cumulative_frozen
    return float(max(0.0, additional_freeze))


def evaluate_partial_melt_step_target(
    row: pd.Series,
    reference_mass: float,
    target_f_melt: float,
    freeze_solid: bool,
    f_melt_tol: float,
) -> tuple[bool, float]:
    """Return whether a target should be skipped plus the current remaining melt fraction."""
    freeze_mass_tol = max(1.0e-9, f_melt_tol * max(reference_mass, 1.0))
    current_remaining_f_melt = compute_remaining_melt_fraction_from_bookkeeping(
        row,
        reference_mass,
    )
    additional_freeze_mass = compute_additional_freeze_mass_from_bookkeeping(
        row,
        reference_mass,
        target_f_melt,
    )
    should_skip = bool(freeze_solid and additional_freeze_mass <= freeze_mass_tol)
    return should_skip, current_remaining_f_melt


def compute_partial_melt_pressure_bar(row: pd.Series) -> float:
    """Reproduce Partial_Melt_Version/calculate_P.cpp in Python for one solved row.

    In the chained partial-melt workflow, ``Pstd`` is used as the current total
    atmosphere pressure carried into the next solve, not as the fixed 1 bar
    standard-state reference used in the standard solver. It therefore needs to
    be recomputed from the evolving atmosphere/interior mass partition every
    step.
    """
    mu = load_atomic_weights()

    moles_atm = float(pd.to_numeric(row.get("Moles_atm", 0.0), errors="coerce"))
    moles_silicate = float(pd.to_numeric(row.get("Moles_silicate", 0.0), errors="coerce"))
    total_active_moles = moles_atm + moles_silicate
    if total_active_moles <= 0.0:
        return float("nan")

    grams_per_mole_atm = 0.0
    for species in (
        "H2_gas",
        "CO_gas",
        "CO2_gas",
        "CH4_gas",
        "O2_gas",
        "H2O_gas",
        "Fe_gas",
        "Mg_gas",
        "SiO_gas",
        "Na_gas",
        "SiH4_gas",
        "SO2_gas",
        "H2S_gas",
        "N2_gas",
        "NH3_gas",
        "HCN_gas",
    ):
        grams_per_mole_atm += float(pd.to_numeric(row.get(species, 0.0), errors="coerce")) * mu.get(species, 0.0)

    grams_per_mole_silicate = 0.0
    for species in SILICATE_COLUMNS:
        grams_per_mole_silicate += float(pd.to_numeric(row.get(species, 0.0), errors="coerce")) * mu.get(species, 0.0)

    mass_atm = (moles_atm / total_active_moles) * grams_per_mole_atm
    mass_silicate = (moles_silicate / total_active_moles) * grams_per_mole_silicate
    frozen_core = float(pd.to_numeric(row.get("M_frozen_core", 0.0), errors="coerce"))
    frozen_solid = float(pd.to_numeric(row.get("M_frozen_solid", 0.0), errors="coerce"))
    total_mass = frozen_core + frozen_solid + mass_atm + mass_silicate
    if total_mass <= 0.0:
        return float("nan")

    matm = mass_atm / total_mass
    interior_fraction = 1.0 - matm
    planet_mass = float(pd.to_numeric(row.get("Planetmass", 0.0), errors="coerce"))
    if interior_fraction <= 0.0 or planet_mass <= 0.0:
        return float("nan")

    return 1.2e6 * matm * np.power(planet_mass * interior_fraction, 2.0 / 3.0) / interior_fraction


def get_GCE_results_state(next_row: dict | pd.Series) -> pd.Series:
    """Return the explicit next-step state used to seed the following solve.

    ``compute_next_partial_melt_state`` stores both the solved previous state and
    the targeted next-step melt inventory in one record. For GCE-state physics such
    as the pressure calculation, we need the targeted next-step state:
    updated active silicate moles/composition plus cumulative frozen mass.
    """
    if isinstance(next_row, pd.Series):
        row = next_row.to_dict()
    else:
        row = dict(next_row)

    gce_state = row.copy()
    gce_state["Moles_silicate"] = float(pd.to_numeric(
        row.get("Moles_silicate_active_next", row.get("Moles_silicate", 0.0)),
        errors="coerce",
    ))
    gce_state["M_frozen_solid"] = float(pd.to_numeric(
        row.get("M_frozen_solid_next", row.get("M_frozen_solid", 0.0)),
        errors="coerce",
    ))

    for species in SILICATE_COLUMNS:
        melt_frac_key = f"{species}_melt_frac_next"
        if melt_frac_key in row:
            gce_state[species] = float(pd.to_numeric(row.get(melt_frac_key, 0.0), errors="coerce"))

    return pd.Series(gce_state)


def build_recorded_partial_melt_step_state(
    solved_row: pd.Series,
    next_step_state: dict[str, float],
    target_f_melt: float,
) -> pd.Series:
    """Return the physically consistent chained state after one solved partial-melt step."""
    current_state = solved_row.copy()
    current_state["f_melt"] = float(target_f_melt)
    current_state["P_SME"] = float(pd.to_numeric(next_step_state.get("P_SME", 0.0), errors="coerce"))
    # Recalculate the varying total atmosphere pressure from the solved state
    # before recording this step. In partial melt, ``Pstd`` evolves with the
    # atmosphere mass fraction rather than staying fixed at 1 bar.
    current_state["Pstd"] = compute_partial_melt_pressure_bar(current_state)

    solved_active_mass = compute_active_silicate_mass(current_state)
    target_active_mass = float(pd.to_numeric(next_step_state.get("M_active_target", 0.0), errors="coerce"))
    current_state["M_active_solved"] = solved_active_mass
    current_state["M_active_target"] = target_active_mass
    current_state["M_active_undershoot"] = max(0.0, target_active_mass - solved_active_mass)
    return current_state


def estimate_minimum_reachable_f_melt(row: pd.Series) -> tuple[float, float, float]:
    """Estimate the lowest reachable f_melt from one state when volatile silicates cannot freeze."""
    mu = load_atomic_weights()
    moles_silicate = float(pd.to_numeric(row.get("Moles_silicate", 0.0), errors="coerce"))
    nonfreezable_mass = 0.0
    for species in VOLATILE_SILICATE_COLUMNS:
        if species not in row.index:
            continue
        species_fraction = float(pd.to_numeric(row.get(species, 0.0), errors="coerce"))
        nonfreezable_mass += species_fraction * moles_silicate * mu.get(species, 0.0)

    m_silicate_ref = compute_active_silicate_mass(row)
    if m_silicate_ref <= 0.0:
        return 0.0, nonfreezable_mass, m_silicate_ref
    return nonfreezable_mass / m_silicate_ref, nonfreezable_mass, m_silicate_ref


# Core-freezing helpers
# These reinterpret a solved full-melt metal phase as a frozen core reservoir.
def _compute_frozen_core_mass(df: pd.DataFrame, mol_weights: dict[str, float]) -> pd.Series:
    """Estimate the frozen core mass from the solved metal phase of each full-melt case."""
    if "Moles_metal" not in df.columns:
        return pd.Series(0.0, index=df.index, dtype=float)

    grams_per_mole_metal = pd.Series(0.0, index=df.index, dtype=float)
    for species in METAL_COLUMNS:
        if species not in df.columns:
            continue
        grams_per_mole_metal = grams_per_mole_metal + df[species].astype(float) * mol_weights[species]

    moles_atm = pd.to_numeric(df["Moles_atm"], errors="coerce").fillna(0.0).astype(float) if "Moles_atm" in df.columns else pd.Series(0.0, index=df.index, dtype=float, name="Moles_atm")
    moles_silicate = pd.to_numeric(df["Moles_silicate"], errors="coerce").fillna(0.0).astype(float) if "Moles_silicate" in df.columns else pd.Series(0.0, index=df.index, dtype=float, name="Moles_silicate")
    moles_metal = pd.to_numeric(df["Moles_metal"], errors="coerce").fillna(0.0).astype(float) if "Moles_metal" in df.columns else pd.Series(0.0, index=df.index, dtype=float, name="Moles_metal")
    moles_total = moles_atm + moles_silicate + moles_metal
    metal_fraction = moles_metal / moles_total.where(moles_total != 0.0, 1.0)
    return metal_fraction * grams_per_mole_metal


def add_frozen_core_columns(df: pd.DataFrame, mol_weights: dict[str, float]) -> pd.DataFrame:
    """Add frozen-core elemental inventories and active pre-split budgets to each row."""
    df = df.copy()
    phase_moles = {
        "metal": pd.to_numeric(df["Moles_metal"], errors="coerce").fillna(0.0).astype(float)
        if "Moles_metal" in df.columns
        else pd.Series(0.0, index=df.index, dtype=float, name="Moles_metal")
    }
    for element in ELEMENTS:
        frozen = pd.Series(0.0, index=df.index, dtype=float)
        for species, coeff in ELEMENT_SPECIES.get(element, []):
            if species not in METAL_COLUMNS or species not in df.columns:
                continue
            frozen = frozen + coeff * df[species].astype(float) * phase_moles["metal"]
        df[f"n{element}_frozen_core"] = frozen

    df["M_frozen_core"] = _compute_frozen_core_mass(df, mol_weights)

    for element in ELEMENTS:
        total_col = f"n{element}"
        frozen_col = f"n{element}_frozen_core"
        active_col = f"n{element}_active_pre_melt_split"
        total = pd.to_numeric(df[total_col], errors="coerce").fillna(0.0).astype(float) if total_col in df.columns else pd.Series(0.0, index=df.index, dtype=float, name=total_col)
        frozen = pd.to_numeric(df[frozen_col], errors="coerce").fillna(0.0).astype(float) if frozen_col in df.columns else pd.Series(0.0, index=df.index, dtype=float, name=frozen_col)
        df[active_col] = total - frozen
    return df


# Step-to-step melt evolution
# This is the main partial-melt state update: given one solved melt state and a
# new target f_melt, compute what remains active in the melt and what freezes out.
def compute_next_partial_melt_state(
    previous_row: pd.Series,
    reference_record: dict[str, float],
    target_f_melt: float,
    volatile_retention_in_solid: bool,
    freeze_solid: bool,
) -> dict[str, float]:
    """Compute the next chained partial-melt state from one solved melt state.

    Responsibilities:
    - decide the target active melt mass for the requested `f_melt`
    - keep volatile silicates in the melt when configured to do so
    - split each silicate species into melt and frozen-solid contributions
    - update cumulative frozen-solid mass and active elemental inventories

    This returns the next state dictionary only; writing `chem_input.dat` and
    saving files happens elsewhere.
    """
    mol_weights = load_atomic_weights()

    next_row = previous_row.to_dict()
    value = pd.to_numeric(previous_row.get("f_melt", 1.0), errors="coerce")
    next_row["f_melt_previous"] = 1.0 if pd.isna(value) else float(value)
    next_row["f_melt_target"] = float(target_f_melt)
    next_row["volatile_retention_in_solid"] = bool(volatile_retention_in_solid)
    next_row["M_silicate_ref"] = float(reference_record["M_silicate_ref"])
    final_full_solidification = freeze_solid and float(target_f_melt) <= 1.0e-12

    species_total_moles: dict[str, float] = {}
    species_masses: dict[str, float] = {}
    current_active_mass = 0.0
    nonfreezable_mass = 0.0
    freezable_mass = 0.0

    # First summarize how much of the current silicate reservoir can freeze at all.
    value = pd.to_numeric(previous_row.get("Moles_silicate", 0.0), errors="coerce")
    current_moles_silicate = 0.0 if pd.isna(value) else float(value)
    for species in SILICATE_COLUMNS:
        value = pd.to_numeric(previous_row.get(species, 0.0), errors="coerce")
        species_total_moles[species] = (0.0 if pd.isna(value) else float(value)) * current_moles_silicate
        species_masses[species] = species_total_moles[species] * mol_weights.get(species, 0.0)
        current_active_mass += species_masses[species]
        if species in VOLATILE_SILICATE_COLUMNS and not volatile_retention_in_solid and not final_full_solidification:
            nonfreezable_mass += species_masses[species]
        else:
            freezable_mass += species_masses[species]

    next_row["M_active_previous"] = current_active_mass

    # Then determine how much additional solid should freeze in this step.
    # The chain schedule is defined against the original reference melt mass, so
    # we freeze enough to reach the desired cumulative frozen mass, not enough
    # to force the post-rebalance silicate phase itself to equal
    # target_f_melt * M_ref.
    if final_full_solidification:
        target_active_mass = 0.0
        mass_to_freeze = current_active_mass
        freeze_fraction = 1.0
    elif not freeze_solid:
        target_active_mass = current_active_mass
        mass_to_freeze = 0.0
        freeze_fraction = 0.0
    else:
        requested_freeze_mass = compute_additional_freeze_mass_from_bookkeeping(
            previous_row,
            float(reference_record["M_silicate_ref"]),
            float(target_f_melt),
        )
        mass_to_freeze = min(max(requested_freeze_mass, 0.0), freezable_mass)
        target_active_mass = current_active_mass - mass_to_freeze
        freeze_fraction = 0.0 if freezable_mass <= 0.0 else mass_to_freeze / freezable_mass

    next_row["M_active_target"] = target_active_mass
    next_row["M_frozen_increment"] = mass_to_freeze

    # Convert that mass-level decision into species-by-species melt and frozen-solid inventories.
    melt_species_moles: dict[str, float] = {}
    solid_species_moles: dict[str, float] = {}
    active_silicate_moles_next = 0.0
    frozen_increment_mass = 0.0

    for species in SILICATE_COLUMNS:
        total_moles = species_total_moles[species]
        if final_full_solidification:
            solid_moles = total_moles
            melt_moles = 0.0
        elif species in VOLATILE_SILICATE_COLUMNS and not volatile_retention_in_solid:
            solid_moles = 0.0
            melt_moles = total_moles
        else:
            solid_moles = freeze_fraction * total_moles
            melt_moles = total_moles - solid_moles
        melt_species_moles[species] = melt_moles
        solid_species_moles[species] = solid_moles
        active_silicate_moles_next += melt_moles
        frozen_increment_mass += solid_moles * mol_weights.get(species, 0.0)

    next_row["Moles_silicate_active_next"] = active_silicate_moles_next
    value = pd.to_numeric(previous_row.get("M_frozen_solid", 0.0), errors="coerce")
    next_row["M_frozen_solid_previous"] = 0.0 if pd.isna(value) else float(value)
    next_row["M_frozen_solid_next"] = next_row["M_frozen_solid_previous"] + frozen_increment_mass

    # Finally update the active elemental inventory after removing what froze out.
    frozen_element_increment = {element: 0.0 for element in ELEMENTS}
    for element in ELEMENTS:
        for species, coeff in ELEMENT_SPECIES.get(element, []):
            if species in solid_species_moles:
                frozen_element_increment[element] += coeff * solid_species_moles[species]

    for element in ELEMENTS:
        value = pd.to_numeric(previous_row.get(f"n{element}", 0.0), errors="coerce")
        current_active = 0.0 if pd.isna(value) else float(value)
        next_row[f"n{element}_frozen_increment"] = frozen_element_increment[element]
        next_row[f"n{element}_active_next"] = current_active - frozen_element_increment[element]

    for species in SILICATE_COLUMNS:
        melt_moles = melt_species_moles[species]
        solid_moles = solid_species_moles[species]
        next_row[f"{species}_melt_moles"] = melt_moles
        next_row[f"{species}_solid_moles"] = solid_moles
        next_row[f"{species}_melt_frac_next"] = (
            melt_moles / active_silicate_moles_next if active_silicate_moles_next > 0.0 else 0.0
        )

    # Store the updated total atmosphere pressure for the next chained solve.
    # The next chem_input.dat should see the pressure implied by the new mass
    # partition after freezing/repartitioning.
    next_row["Pstd"] = compute_partial_melt_pressure_bar(get_GCE_results_state(next_row))

    return next_row
