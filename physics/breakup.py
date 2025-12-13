# breakup.py
"""
Structural breakup model for meteoroids.

This module implements progressive fragmentation based on:
- Aerodynamic loading (dynamic pressure, bending moments)
- Material strength (Weibull size scaling, phase-dependent)
- Thermal weakening (temperature gradients, phase changes)

The fragmentation process is progressive, not binary:
- Below critical stress: no fragmentation
- Near critical stress: slow fragment shedding
- Above critical stress: rapid cascading breakup
"""
import math
import random
from dataclasses import dataclass, field


@dataclass
class StrengthParams:
    """Material strength parameters."""
    base_tensile: float = 2.0e7   # Pa (reference tensile strength)
    weibull_m: float = 7.0        # Weibull shape parameter (higher = less scatter)
    ref_radius: float = 0.5       # m (reference size for Weibull scaling)
    # Phase-dependent strength multiplier
    phase_mult: dict = field(default_factory=lambda: {
        "solid": 1.0,
        "fusion": 0.6,
        "liquid": 0.3,
        "sublimation": 0.5,
        "vapor": 0.1
    })


@dataclass
class AeroParams:
    """Aerodynamic loading parameters."""
    bend_coeff: float = 0.8       # Bending stress efficiency
    shape_factor: float = 1.0     # Bluntness correction (1.0 for sphere)


@dataclass
class FragmentationState:
    """Tracks progressive fragmentation state."""
    stress_integral: float = 0.0   # Accumulated stress-time
    fragments_shed: int = 0        # Total fragments shed so far
    last_stress_ratio: float = 0.0 # Most recent stress ratio


class BreakupModel:
    def __init__(self, strength: StrengthParams | None = None,
                 aero: AeroParams | None = None):
        self.S = strength or StrengthParams()
        self.A = aero or AeroParams()
        self.state = FragmentationState()

        # Fragmentation thresholds
        self.onset_ratio = 0.7      # Start shedding fragments above this
        self.cascade_ratio = 1.5    # Rapid breakup above this
        self.max_fragment_rate = 20 # Max fragments per timestep

    def effective_tensile(self, props: dict, radius: float, phase: str,
                           thermal_grad_K: float) -> float:
        """
        Calculate effective tensile strength with all degradation factors.

        Uses Weibull scaling: larger bodies have lower effective strength
        due to increased probability of critical flaws.

        Parameters
        ----------
        props : dict
            Material properties
        radius : float
            Body radius [m]
        phase : str
            Current phase
        thermal_grad_K : float
            Surface-core temperature difference [K]

        Returns
        -------
        float
            Effective tensile strength [Pa]
        """
        base = props.get("tensile", self.S.base_tensile)

        # Weibull size scaling: σ_eff = σ_0 * (r_0/r)^(1/m)
        size_mult = (self.S.ref_radius / max(1e-3, radius)) ** (1.0 / max(1.0, self.S.weibull_m))

        # Phase-dependent weakening
        phase_mult = self.S.phase_mult.get(phase, 1.0)

        # Thermal stress weakening
        # Linear penalty: 25% reduction per 1000 K gradient
        thermal_penalty = max(0.3, 1.0 - 0.25 * (thermal_grad_K / 1000.0))

        return max(1e5, base * size_mult * phase_mult * thermal_penalty)

    def aerodynamic_stress(self, rho: float, v: float, radius: float,
                            angle_rad: float) -> tuple:
        """
        Calculate aerodynamic stresses on the body.

        Returns both dynamic pressure and bending stress components.

        Parameters
        ----------
        rho : float
            Atmospheric density [kg/m³]
        v : float
            Velocity [m/s]
        radius : float
            Body radius [m]
        angle_rad : float
            Entry angle from vertical [rad]

        Returns
        -------
        tuple
            (q_dynamic, sigma_bending) in Pa
        """
        # Dynamic pressure
        q_dyn = 0.5 * rho * v * v

        # Bending stress from asymmetric loading
        # Scales with attack angle and body size
        aoa_factor = max(0.2, abs(math.sin(angle_rad)))
        size_factor = radius / max(1e-3, self.S.ref_radius)
        sigma_bend = self.A.bend_coeff * q_dyn * aoa_factor * self.A.shape_factor * size_factor

        return q_dyn, sigma_bend

    def stress_ratio(self, meteor, rho: float, v: float, angle_rad: float) -> float:
        """
        Calculate the stress ratio (applied stress / strength).

        Values > 1.0 indicate structural failure.

        Returns
        -------
        float
            Stress ratio (dimensionless)
        """
        if meteor.is_fully_ablated():
            return 0.0

        thermal_grad = max(0.0, meteor.temperature - meteor.core_temperature)
        sigma_eff = self.effective_tensile(
            meteor.props, max(1e-6, meteor.radius), meteor.phase, thermal_grad
        )

        q_dyn, sigma_bend = self.aerodynamic_stress(
            rho, v, max(1e-6, meteor.radius), angle_rad
        )

        # Combined stress (conservative addition)
        sigma_total = sigma_bend + q_dyn
        ratio = sigma_total / max(1e-6, sigma_eff)

        self.state.last_stress_ratio = ratio
        return ratio

    def should_breakup(self, meteor, rho: float, v: float, angle_rad: float) -> bool:
        """
        Determine if fragmentation should occur this timestep.

        Uses progressive fragmentation: probability increases with stress ratio.
        """
        if meteor.is_fully_ablated():
            return False

        ratio = self.stress_ratio(meteor, rho, v, angle_rad)

        # Below onset: no fragmentation
        if ratio < self.onset_ratio:
            return False

        # Above cascade: always fragment
        if ratio > self.cascade_ratio:
            return True

        # In transition zone: probabilistic based on stress excess
        excess = (ratio - self.onset_ratio) / (self.cascade_ratio - self.onset_ratio)
        probability = excess ** 2  # Quadratic probability curve
        return random.random() < probability

    def compute_fragment_count(self, meteor, rho: float, v: float,
                                 angle_rad: float) -> int:
        """
        Determine how many fragments to shed this timestep.

        Fragment count scales with:
        - Stress excess above threshold
        - Body mass (larger bodies fragment into more pieces)

        Returns
        -------
        int
            Number of fragments to create (0 if no fragmentation)
        """
        ratio = self.state.last_stress_ratio

        if ratio < self.onset_ratio:
            return 0

        # Normalized stress excess
        if ratio > self.cascade_ratio:
            excess = 1.0
        else:
            excess = (ratio - self.onset_ratio) / (self.cascade_ratio - self.onset_ratio)

        # Base fragment count scales with mass (log scale)
        mass_factor = math.log10(max(1.0, meteor.mass)) + 1.0  # 1-4 for 1kg-1000kg

        # More fragments at higher stress
        base_count = 1 + int(excess * mass_factor * 3)

        # Add some randomness
        count = base_count + random.randint(0, max(1, base_count // 2))

        return min(count, self.max_fragment_rate)

    def fragment_cascade(self, meteor, n_fragments: int = None):
        """
        Create fragments from the meteor body.

        Uses power-law mass distribution (smaller fragments more common).

        Parameters
        ----------
        meteor : Meteor
            The meteor object to fragment
        n_fragments : int, optional
            Number of fragments to create. If None, uses compute_fragment_count.
        """
        if meteor.mass <= 0.0:
            return

        if n_fragments is None:
            n_fragments = random.randint(4, 16)

        # Limit total fragments to reasonable number
        existing = len(meteor.particles) if meteor.particles else 0
        max_total = 100
        n_fragments = min(n_fragments, max_total - existing)

        if n_fragments <= 0:
            return

        # Mass to convert to fragments (fraction of current mass)
        ratio = self.state.last_stress_ratio
        if ratio >= self.cascade_ratio:
            fragment_fraction = 0.5 + 0.4 * random.random()  # 50-90% in cascade
        else:
            fragment_fraction = 0.1 + 0.2 * random.random()  # 10-30% normally

        fragment_mass_total = meteor.mass * fragment_fraction
        masses = self._sample_powerlaw_masses(fragment_mass_total, n_fragments, alpha=1.8)

        # Initialize fragmentation if first time
        if not meteor.fragmented:
            meteor.note_fragmentation()

        # Add new particles
        if meteor.particles:
            density = meteor.props["density"]
            for i, m in enumerate(masses):
                if m < 1e-6:  # Skip tiny fragments
                    continue

                # Position fragments trailing the core
                xi = random.random()
                x_pos = -meteor.radius * (1.0 + 2.0 * xi)

                from entities.particles import Particle
                p = Particle(
                    x=x_pos,
                    r=0.01,  # Will be set by shrink_from_mass
                    m=m,
                    T=0.7 * meteor.temperature + 0.3 * meteor.core_temperature,
                    Cd=meteor.props["Cd"] * (1.0 + 0.2 * random.random()),
                    emissivity=meteor.props["emissivity"],
                    cp=meteor.props["cp"],
                    L=meteor.props.get("L_vapor", meteor.props.get("L", 1e6)),
                )
                p.shrink_from_mass(density)
                p.phase = meteor.phase
                meteor.particles.append(p)

        # Reduce core mass
        meteor.mass = max(0.0, meteor.mass - fragment_mass_total)
        meteor.update_geometry_from_mass()

        # Update tracking
        self.state.fragments_shed += n_fragments

    def _sample_powerlaw_masses(self, total_mass: float, N: int, alpha: float = 1.8) -> list:
        """
        Sample N masses from a power-law distribution.

        Uses inverse transform sampling on m ~ x^(-alpha).
        Normalized to sum to total_mass.
        """
        if N <= 0 or total_mass <= 0:
            return []

        x_min, x_max = 1.0, 100.0
        samples = []

        for _ in range(N):
            u = random.random()
            x = ((u * (x_max**(1.0 - alpha) - x_min**(1.0 - alpha)) +
                  x_min**(1.0 - alpha)) ** (1.0 / (1.0 - alpha)))
            samples.append(x)

        total = sum(samples)
        if total < 1e-9:
            return [total_mass / N] * N

        return [total_mass * (x / total) for x in samples]

    def reset_state(self):
        """Reset fragmentation tracking state."""
        self.state = FragmentationState()
