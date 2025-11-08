# vaporization_detector.py
from collections import deque

class VaporizationDetector:
    def __init__(self, window=25):  # ~0.5s if dt=0.02
        self.window = window
        self.history = deque(maxlen=window)

    def update(self, phase: str):
        self.history.append(phase == "vapor")

    def sustained(self, threshold=0.6) -> bool:
        if not self.history:
            return False
        frac = sum(self.history) / len(self.history)
        return frac >= threshold