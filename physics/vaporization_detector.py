# vaporization_detector.py
"""
Detector for sustained vaporization/sublimation events.

This module tracks the phase and temperature history to determine if the material
has been in a destructive (vaporizing) state for a sustained period at HIGH
temperatures. This indicates actual core destruction, not just thermodynamic
tendency to sublimate.

Key distinction:
- Thermodynamic "sublimation" phase at low T: negligible mass loss
- Sustained vapor phase at high T (near fusion): actual destruction
"""
from collections import deque


class VaporizationDetector:
    def __init__(self, window: int = 50):
        """
        Initialize the detector.

        Parameters
        ----------
        window : int
            Number of timesteps to track (~1s if dt=0.02)
        """
        self.window = window
        self.phase_history = deque(maxlen=window)
        self.temp_history = deque(maxlen=window)
        self.fusion_T = 1000.0  # Will be updated

    def update(self, phase: str, T_surface: float, fusion_T: float):
        """
        Record the current phase and temperature.

        Parameters
        ----------
        phase : str
            Current phase: "solid", "fusion", "liquid", "sublimation", or "vapor"
        T_surface : float
            Current surface temperature [K]
        fusion_T : float
            Material's fusion temperature [K]
        """
        is_vaporizing = phase in ("vapor", "sublimation")
        self.phase_history.append(is_vaporizing)
        self.temp_history.append(T_surface)
        self.fusion_T = fusion_T

    def sustained(self, threshold: float = 0.6) -> bool:
        """
        Check if vaporization/sublimation has been sustained.
        (Legacy method - just checks phase, not temperature)

        Parameters
        ----------
        threshold : float
            Fraction of recent timesteps that must be in destructive phase

        Returns
        -------
        bool
            True if sustained vaporization/sublimation detected
        """
        if not self.phase_history:
            return False
        frac = sum(self.phase_history) / len(self.phase_history)
        return frac >= threshold

    def sustained_destruction(self, phase_threshold: float = 0.8,
                               temp_threshold_fraction: float = 1.5) -> bool:
        """
        Check if the material is being actively destroyed.

        Requires BOTH:
        1. Sustained vaporizing phase (>80% of recent history)
        2. High temperature (>150% of fusion temperature on average)

        The temperature threshold is set well ABOVE the fusion point to distinguish:
        - Slow sublimation at moderate temperatures (normal, not catastrophic)
        - Rapid vaporization at extreme temperatures (actual destruction)

        For ice (fusion_T=273K), this requires T > 410K (137°C)
        For stone (fusion_T=1700K), this requires T > 2550K

        Parameters
        ----------
        phase_threshold : float
            Fraction of timesteps that must be in vapor/sublimation phase
        temp_threshold_fraction : float
            Average temperature must be this fraction of fusion_T

        Returns
        -------
        bool
            True if core is being destroyed by vaporization
        """
        if len(self.phase_history) < self.window // 2:
            # Not enough history yet
            return False

        # Check phase criterion
        phase_frac = sum(self.phase_history) / len(self.phase_history)
        if phase_frac < phase_threshold:
            return False

        # Check temperature criterion
        # Must be significantly above fusion point to indicate destruction
        avg_temp = sum(self.temp_history) / len(self.temp_history)
        temp_threshold = temp_threshold_fraction * self.fusion_T

        return avg_temp >= temp_threshold

    def reset(self):
        """Clear the history."""
        self.phase_history.clear()
        self.temp_history.clear()
