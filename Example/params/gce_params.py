from dataclasses import dataclass, field
from typing import Dict
import numpy as np


@dataclass
class GCEParams:
    Planetmassarray: np.ndarray = field(default_factory=lambda: np.array([5.]))
    FakeMolesTotal: float = 10e3
    T_AMOI_array: np.ndarray = field(default_factory=lambda: np.array([2500.]))
    T_SME_array: np.ndarray = field(default_factory=lambda: np.array([3000.]))
    tarmgsiarray: np.ndarray = field(default_factory=lambda: np.array([1.]))
    tarfesiarray: np.ndarray = field(default_factory=lambda: np.array([1.]))
    tarWaterarray: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    tarHHearray: np.ndarray = field(default_factory=lambda: np.array([0.03]))
    tarDiskCOarray: np.ndarray = field(default_factory=lambda: np.array([0.5]))
    tarDiskSHarray: np.ndarray = field(default_factory=lambda: np.array([1.335e-5]))
    tarFearray: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    tarOarray: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    P_AMOI_array: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    P_SME_array: np.ndarray = field(default_factory=lambda: np.array([10.0]))
    UseCondriticComp: str = "molar fraction"
    UseCondriticPreset: str = "ed_young"


    condritic_mass_allegre: Dict[str, float] = field(
        default_factory=lambda: {
            "Si": 0.171,
            "Mg": 0.158,
            "Fe": 0.288,
            "Na": 0.00187,
            "O": 0.32436,
            "C": 0.0017,
            "S": 0.0,
        }
    )
    condritic_mass_javoy: Dict[str, float] = field(
        default_factory=lambda: {
            "Si": 0.1923,
            "Mg": 0.1221,
            "Fe": 0.3339,
            "Na": 0.00187,
            "O": 0.3028,
            "C": 0.0010,
            "S": 0.0,
        }
    )
    condritic_molar_ed_young: Dict[str, float] = field( # from Kallemeyn & Wasson 1986 and Grady et al. 1986
        default_factory=lambda: {
            "O": 0.4962,
            "Mg": 0.167,
            "Si": 0.164,
            "C": 0.011,
            "Fe": 0.159,
            "Na": 0.0028,
            # almost matches except for O and Fe
            "S": 0.0236, # from Lodders 2021: this is sadly CI chondrite
            "N": 0.00239, # from Lodders 2021: this is sadly CI chondrite
        }
    )

    def select_condritic(self, sulfur_enabled: bool, nitrogen_enabled: bool = False) -> Dict[str, float]:
        """Select and return chondritic composition, zeroing S/N if their versions are disabled."""
        if self.UseCondriticComp == "mass fraction":
            if self.UseCondriticPreset == "allegre":
                cond = self.condritic_mass_allegre
            elif self.UseCondriticPreset == "javoy":
                cond = self.condritic_mass_javoy
            else:
                raise ValueError(
                    "For mass fraction, UseCondriticPreset must be 'allegre' or 'javoy'"
                )
        elif self.UseCondriticComp == "molar fraction":
            cond = self.condritic_molar_ed_young
        else:
            raise ValueError("UseCondriticComp must be 'mass fraction' or 'molar fraction'")

        if not sulfur_enabled or not nitrogen_enabled:
            cond = cond.copy()
            if not sulfur_enabled:
                cond["S"] = 0.0
            if not nitrogen_enabled:
                cond["N"] = 0.0

        return cond
