# vaporization_detector.py
"""
Detector for sustained vaporization/sublimation events.

This module tracks the phase history and determines if the material
has been in a vaporizing or sublimating state for a sustained period,
indicating potential complete destruction of the meteoroid.
"""
from collections import deque


class VaporizationDetector:
    def __init__(self, window: int = 25):
        """
        Initialize the detector.

        Parameters
        ----------
        window : int
            Number of timesteps to track (~0.5s if dt=0.02)
        """
        self.window = window
        self.history = deque(maxlen=window)

    def update(self, phase: str):
        """
        Record the current phase.

        Parameters
        ----------
        phase : str
            Current phase: "solid", "fusion", "liquid", "sublimation", or "vapor"
        """
        # Track both vaporization and sublimation as "destructive" phases
        is_destructive = phase in ("vapor", "sublimation")
        self.history.append(is_destructive)

    def sustained(self, threshold: float = 0.6) -> bool:
        """
        Check if vaporization/sublimation has been sustained.

        Parameters
        ----------
        threshold : float
            Fraction of recent timesteps that must be in destructive phase

        Returns
        -------
        bool
            True if sustained vaporization/sublimation detected
        """
        if not self.history:
            return False
        frac = sum(self.history) / len(self.history)
        return frac >= threshold

    def reset(self):
        """Clear the phase history."""
        self.history.clear()