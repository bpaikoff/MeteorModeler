# physics/__init__.py
"""Physics models for meteor atmospheric entry simulation."""

from physics.heating import HeatingModel, HeatingCoeffs, ThermalProperties, THERMAL_CONDUCTIVITY
from physics.ablation import AblationModel
from physics.phase import PhaseModel
from physics.vapor_pressure import VaporPressureModel
from physics.energy_balance import EnergyBalanceModel, EnergyPartition
from physics.breakup import BreakupModel, StrengthParams, AeroParams, FragmentationState
from physics.vaporization_detector import VaporizationDetector
from physics.atmosphere import Atmosphere
from physics.trajectory import Trajectory

__all__ = [
    "HeatingModel", "HeatingCoeffs", "ThermalProperties", "THERMAL_CONDUCTIVITY",
    "AblationModel",
    "PhaseModel",
    "VaporPressureModel",
    "EnergyBalanceModel", "EnergyPartition",
    "BreakupModel", "StrengthParams", "AeroParams", "FragmentationState",
    "VaporizationDetector",
    "Atmosphere",
    "Trajectory",
]
