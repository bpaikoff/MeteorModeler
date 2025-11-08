# heating.py
import math
from dataclasses import dataclass

@dataclass
class HeatingCoeffs:
    # Sutton-Graves-like convective coefficient (simplified units)
    k_sg: float = 1.83e-4
    # Fay-Riddell augmentation factor
    k_fr: float = 3.5e-4
    # Plasma/shock-layer radiation coefficient
    k_rad: float = 1.2e-4
    # Velocity exponents (radiation steeper)
    b_conv: float = 3.0
    b_rad: float = 2.9
    # Density exponents
    a_conv: float = 0.5
    a_rad: float = 0.8
    # Catalytic recombination multiplier (0–1)
    catalytic: float = 0.2
    # Shape/bluntness factor (caps heating at stagnation)
    shape_factor: float = 1.0

class HeatingModel:
    def __init__(self, coeffs: HeatingCoeffs | None = None):
        self.sigma = 5.670374419e-8
        self.coeffs = coeffs or HeatingCoeffs()
        # Regime threshold by altitude; switch blends around ~90–100 km
        self.fm_threshold_alt_m = 95_000.0
        self.fm_halfwidth_m = 15_000.0  # smooth blend range

    def _blend_free_molecular(self, altitude_m: float) -> float:
        # Smooth logistic blend: 1 at high alt, 0 in continuum
        x = (altitude_m - self.fm_threshold_alt_m) / self.fm_halfwidth_m
        return 1.0 / (1.0 + math.exp(-x))

    def emissivity_T(self, base_eps: float, T_surface: float, composition: str) -> float:
        # Temperature-dependent emissivity; carbonaceous/stone radiate efficiently at high T
        if composition in ("carbonaceous", "stone"):
            return min(0.98, base_eps + 0.05 * (min(T_surface, 3000.0) - 1000.0) / 2000.0)
        # Ice/iron stay closer to base with mild increase
        return min(0.95, base_eps + 0.03 * (min(T_surface, 2500.0) - 800.0) / 1700.0)

    def q_convective(self, rho: float, v: float, blend_fm: float) -> float:
        # Blend free-molecular and continuum convective heating
        k_fr = self.coeffs.k_fr
        k_sg = self.coeffs.k_sg
        a = self.coeffs.a_conv
        b = self.coeffs.b_conv

        q_cont = k_sg * (rho ** a) * (v ** b)
        q_fm = 0.6 * k_fr * (rho ** 0.3) * (v ** (b - 0.2))  # FM weaker rho dependence
        return (1.0 - blend_fm) * q_cont + blend_fm * q_fm

    def q_plasma(self, rho: float, v: float, blend_fm: float) -> float:
        # Shock-layer radiation; reduced in FM regime
        k = self.coeffs.k_rad
        a = self.coeffs.a_rad
        b = self.coeffs.b_rad
        q_cont = k * (rho ** a) * (v ** b)
        q_fm = 0.5 * k * (rho ** 0.4) * (v ** (b - 0.2))
        return (1.0 - blend_fm) * q_cont + blend_fm * q_fm

    def q_surface_radiative_cooling(self, emissivity: float, T_surface: float, T_air: float) -> float:
        return emissivity * self.sigma * (T_surface**4 - T_air**4)

    def net_heat_flux(self, altitude_m: float, rho: float, v: float,
                       emissivity_base: float, T_surface: float, T_air: float,
                       composition: str, angle_rad: float, fragment_view_factor: float = 0.0) -> float:
        # Regime blending
        blend_fm = self._blend_free_molecular(altitude_m)
        # Angle-of-attack diminishes stagnation flux
        aoa_factor = max(0.3, math.cos(max(0.0, angle_rad)))  # blunt entry gets most heating
        shape = self.coeffs.shape_factor

        q_conv = self.q_convective(rho, v, blend_fm)
        q_plasma = self.q_plasma(rho, v, blend_fm)
        # Catalytic recombination pushes convective term up in oxidizing flows
        q_conv *= (1.0 + self.coeffs.catalytic)

        eps_T = self.emissivity_T(emissivity_base, T_surface, composition)
        q_cool = self.q_surface_radiative_cooling(eps_T, T_surface, T_air)

        # Fragment view factor increases re-radiation onto core (optional)
        q_fragment_reradiation = fragment_view_factor * 0.2 * q_plasma

        q_net = (q_conv + q_plasma + q_fragment_reradiation) * aoa_factor * shape - q_cool
        return max(0.0, q_net)

    def update_surface_and_core_temps(self, altitude_m: float, A: float, m: float, cp: float,
                                      emissivity: float, rho: float, v: float,
                                      T_surface: float, T_core: float, T_air: float,
                                      dt: float, composition: str, angle_rad: float,
                                      fragment_view_factor: float = 0.0, core_conduction_fraction: float = 0.10):
        q_net = self.net_heat_flux(altitude_m, rho, v, emissivity, T_surface, T_air,
                                   composition, angle_rad, fragment_view_factor)
        P_net = q_net * max(1e-9, A)
        dTdt_surface = P_net / (max(1e-9, m) * max(1e-9, cp))

        T_surface += dTdt_surface * dt
        T_core += core_conduction_fraction * dTdt_surface * dt
        return T_surface, T_core, P_net

    def ablation_rate(self, P_net: float, L_vap: float, T_surface: float, fusion_T: float, vapor_T: float) -> float:
        if T_surface < 0.9 * fusion_T:
            return 0.0
        superheat = max(0.0, T_surface - fusion_T)
        if T_surface < vapor_T:
            f = 0.15 + 0.25 * (superheat / (superheat + 500.0))
        else:
            f = 0.5 + 0.35 * (superheat / (superheat + 1000.0))
        f = min(0.9, f)
        return (f * max(0.0, P_net)) / max(1e-9, L_vap)