import math
from dataclasses import dataclass

@dataclass
class Particle:
    # Minimal state per fragment/parcel
    x: float        # along-path coordinate (m); 0 = leading edge, positive trailing
    r: float        # effective radius (m)
    m: float        # mass (kg)
    T: float        # temperature (K)
    Cd: float       # drag coefficient (per fragment)
    emissivity: float
    cp: float       # J/(kg*K)
    L: float        # J/kg (ablation energy)
    alive: bool = True

    @property
    def A(self) -> float:
        return math.pi * self.r * self.r

    def shrink_from_mass(self, density: float):
        # update radius from mass assuming sphere
        if self.m <= 0:
            self.r = 0.0
            self.alive = False
            return
        V = self.m / density
        self.r = ((3.0 * V) / (4.0 * math.pi)) ** (1.0 / 3.0)