# breakup.py
import math
import random
from dataclasses import dataclass

@dataclass
class StrengthParams:
    base_tensile: float = 2.0e7   # Pa
    weibull_m: float = 7.0        # shape parameter
    ref_radius: float = 0.5       # m
    phase_mult = {"solid": 1.0, "fusion": 0.6, "liquid": 0.3, "vapor": 0.1}

@dataclass
class AeroParams:
    bend_coeff: float = 0.8       # aerodynamic bending efficiency
    shape_factor: float = 1.0     # bluntness correction

class BreakupModel:
    def __init__(self, strength: StrengthParams | None = None, aero: AeroParams | None = None):
        self.S = strength or StrengthParams()
        self.A = aero or AeroParams()

    def effective_tensile(self, props: dict, radius: float, phase: str, thermal_grad_K: float) -> float:
        # Weibull scaling: larger bodies have lower effective strength
        base = props.get("tensile", self.S.base_tensile)
        size_mult = (self.S.ref_radius / max(1e-3, radius)) ** (1.0 / max(1.0, self.S.weibull_m))
        phase_mult = self.S.phase_mult.get(phase, 1.0)
        # Thermal stress reduces strength; linear penalty per 1000 K gradient
        thermal_penalty = max(0.5, 1.0 - 0.25 * (thermal_grad_K / 1000.0))
        return max(1e5, base * size_mult * phase_mult * thermal_penalty)

    def aerodynamic_bending_stress(self, rho: float, v: float, radius: float, angle_rad: float) -> float:
        # Approximate bending stress from pressure differential across body
        q_dyn = 0.5 * rho * v * v
        # Lever arm ~ radius, correction by angle-of-attack and shape
        aoa = max(0.2, math.sin(max(0.0, angle_rad)))
        sigma_bend = self.A.bend_coeff * q_dyn * aoa * self.A.shape_factor
        # Scale with size (larger bodies experience larger bending moments)
        return sigma_bend * (radius / self.S.ref_radius)

    def should_breakup(self, meteor, rho: float, v: float, angle_rad: float) -> bool:
        if meteor.is_fully_ablated():
            return False
        thermal_grad = max(0.0, meteor.temperature - meteor.core_temperature)
        sigma_eff = self.effective_tensile(meteor.props, max(1e-6, meteor.radius), meteor.phase, thermal_grad)
        sigma_bend = self.aerodynamic_bending_stress(rho, v, max(1e-6, meteor.radius), angle_rad)

        # Additionally, consider dynamic pressure vs strength
        q_dyn = 0.5 * rho * v * v
        stress_ratio = (sigma_bend + q_dyn) / max(1e-6, sigma_eff)
        return stress_ratio > 1.0

    def fragment_cascade(self, meteor):
        # Create a distribution of child fragments; conserve mass approximately
        total_mass = meteor.mass
        if total_mass <= 0.0:
            return
        N = random.randint(32, 96)
        masses = self._sample_powerlaw_masses(total_mass, N, alpha=1.8)
        # Use meteor.note_fragmentation to spawn particles; replace masses uniformly
        meteor.note_fragmentation()
        if meteor.particles:
            alive = [p for p in meteor.particles if p.alive]
            for i, p in enumerate(alive[:len(masses)]):
                p.m = masses[i]
                p.shrink_from_mass(meteor.props["density"])
                # Initialize fragment temperature near outer temperature
                p.T = 0.8 * meteor.temperature + 0.2 * meteor.core_temperature
                p.phase = meteor.phase

    def _sample_powerlaw_masses(self, total_mass: float, N: int, alpha: float = 2.0):
        # Sample N masses ~ power-law, normalized to total_mass
        # m ~ x^(-alpha); generate via inverse transform on [x_min,x_max]
        x_min, x_max = 1.0, 100.0
        samples = []
        for _ in range(N):
            u = random.random()
            x = ((u * (x_max**(1.0 - alpha) - x_min**(1.0 - alpha)) + x_min**(1.0 - alpha)) ** (1.0 / (1.0 - alpha)))
            samples.append(x)
        s = sum(samples)
        return [total_mass * (x / s) for x in samples]