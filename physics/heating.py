# heating.py
"""
Aerodynamic heating model for meteor entry.

This module implements heating from:
1. Convective heating (Sutton-Graves and Fay-Riddell correlations)
2. Shock-layer radiative heating
3. Surface radiative cooling (Stefan-Boltzmann)
4. Internal heat conduction (thermal diffusivity-based)

The model transitions between free-molecular and continuum flow regimes
based on altitude (Knudsen number proxy).
"""
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from physics.plasma import HeatFluxBreakdown


@dataclass
class HeatingCoeffs:
    """Empirical coefficients for heating correlations."""
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


@dataclass
class ThermalProperties:
    """Thermal properties for heat conduction calculations."""
    # Default values for typical stony material
    thermal_conductivity: float = 2.0  # W/(m·K)
    density: float = 3500.0  # kg/m³
    specific_heat: float = 800.0  # J/(kg·K)

    @property
    def thermal_diffusivity(self) -> float:
        """Thermal diffusivity α = k/(ρ·cp) [m²/s]"""
        return self.thermal_conductivity / (self.density * self.specific_heat)


# Thermal conductivities for different materials [W/(m·K)]
THERMAL_CONDUCTIVITY = {
    "iron": 80.0,       # Pure iron ~80 W/(m·K)
    "stone": 2.0,       # Silicates ~1-3 W/(m·K)
    "ice": 2.2,         # Ice ~2.2 W/(m·K)
    "carbonaceous": 1.0 # Porous, low conductivity
}


class HeatingModel:
    def __init__(self, coeffs: HeatingCoeffs | None = None):
        self.sigma = 5.670374419e-8  # Stefan-Boltzmann constant [W/(m²·K⁴)]
        self.coeffs = coeffs or HeatingCoeffs()
        # Regime threshold by altitude; switch blends around ~90–100 km
        self.fm_threshold_alt_m = 95_000.0
        self.fm_halfwidth_m = 15_000.0  # smooth blend range

    def _blend_free_molecular(self, altitude_m: float) -> float:
        """
        Smooth logistic blend between flow regimes.

        Returns 1.0 at high altitude (free-molecular), 0.0 in continuum.
        The transition occurs around 95 km with a 15 km halfwidth.
        """
        x = (altitude_m - self.fm_threshold_alt_m) / self.fm_halfwidth_m
        return 1.0 / (1.0 + math.exp(-x))

    def emissivity_T(self, base_eps: float, T_surface: float, composition: str) -> float:
        """
        Temperature-dependent emissivity.

        Emissivity generally increases with temperature due to increased
        absorption bands and surface oxidation.
        """
        if composition in ("carbonaceous", "stone"):
            # Carbonaceous and stone radiate efficiently at high T
            return min(0.98, base_eps + 0.05 * (min(T_surface, 3000.0) - 1000.0) / 2000.0)
        # Ice/iron stay closer to base with mild increase
        return min(0.95, base_eps + 0.03 * (min(T_surface, 2500.0) - 800.0) / 1700.0)

    def q_convective(self, rho: float, v: float, blend_fm: float) -> float:
        """
        Convective heat flux [W/m²].

        Blends between continuum (Sutton-Graves-like) and free-molecular
        heating based on altitude. Continuum heating scales as:
            q ~ ρ^0.5 * v^3

        Free-molecular heating has weaker density dependence.
        """
        k_fr = self.coeffs.k_fr
        k_sg = self.coeffs.k_sg
        a = self.coeffs.a_conv
        b = self.coeffs.b_conv

        q_cont = k_sg * (rho ** a) * (v ** b)
        q_fm = 0.6 * k_fr * (rho ** 0.3) * (v ** (b - 0.2))
        return (1.0 - blend_fm) * q_cont + blend_fm * q_fm

    def q_plasma(self, rho: float, v: float, blend_fm: float) -> float:
        """
        Shock-layer radiative heat flux [W/m²].

        At high velocities, the shock-heated gas radiates intensely.
        This scales more steeply with density than convective heating.
        Reduced in free-molecular regime where no shock forms.
        """
        k = self.coeffs.k_rad
        a = self.coeffs.a_rad
        b = self.coeffs.b_rad
        q_cont = k * (rho ** a) * (v ** b)
        q_fm = 0.5 * k * (rho ** 0.4) * (v ** (b - 0.2))
        return (1.0 - blend_fm) * q_cont + blend_fm * q_fm

    def q_surface_radiative_cooling(self, emissivity: float, T_surface: float, T_air: float) -> float:
        """
        Surface radiative cooling [W/m²] via Stefan-Boltzmann law.

        q_rad = ε·σ·(T_surface⁴ - T_ambient⁴)

        This is the primary cooling mechanism for meteors at high temperatures.
        """
        return emissivity * self.sigma * (T_surface**4 - T_air**4)

    def equilibrium_temperature(self, q_incident: float, emissivity: float, T_air: float = 0.0) -> float:
        """
        Calculate radiative equilibrium temperature.

        At equilibrium: q_incident = ε·σ·(T⁴ - T_air⁴)
        Solving: T = ((q_incident / (ε·σ)) + T_air⁴)^0.25

        Parameters
        ----------
        q_incident : float
            Incident heat flux [W/m²]
        emissivity : float
            Surface emissivity (0-1)
        T_air : float
            Ambient temperature [K]

        Returns
        -------
        float
            Equilibrium surface temperature [K]
        """
        if q_incident <= 0:
            return T_air

        eps_sigma = max(1e-9, emissivity * self.sigma)
        T4 = q_incident / eps_sigma + T_air**4
        return T4 ** 0.25

    def net_heat_flux(self, altitude_m: float, rho: float, v: float,
                       emissivity_base: float, T_surface: float, T_air: float,
                       composition: str, angle_rad: float, fragment_view_factor: float = 0.0) -> float:
        """
        Net heat flux into the meteor surface [W/m²].

        Delegates to compute_heat_flux_breakdown and returns just the net value.
        Can be negative when radiative cooling exceeds heating.
        """
        breakdown = self.compute_heat_flux_breakdown(
            altitude_m, rho, v, emissivity_base, T_surface, T_air,
            composition, angle_rad, fragment_view_factor
        )
        return breakdown.q_net

    def compute_heat_flux_breakdown(self, altitude_m: float, rho: float, v: float,
                                     emissivity_base: float, T_surface: float, T_air: float,
                                     composition: str, angle_rad: float,
                                     fragment_view_factor: float = 0.0) -> "HeatFluxBreakdown":
        """
        Compute detailed breakdown of heat flux contributions.

        Returns a HeatFluxBreakdown dataclass with individual flux components
        and their relative percentages.
        """
        from physics.plasma import HeatFluxBreakdown

        # Regime blending
        blend_fm = self._blend_free_molecular(altitude_m)

        # Angle-of-attack factor
        aoa_factor = max(0.3, math.cos(max(0.0, angle_rad)))
        shape = self.coeffs.shape_factor

        # Individual flux components
        q_conv_raw = self.q_convective(rho, v, blend_fm)
        q_plasma_raw = self.q_plasma(rho, v, blend_fm)

        # Catalytic augmentation (stored separately for display)
        q_catalytic = q_conv_raw * self.coeffs.catalytic

        # Total convective with catalytic
        q_conv_total = q_conv_raw * (1.0 + self.coeffs.catalytic)

        # Fragment re-radiation
        q_fragment = fragment_view_factor * 0.2 * q_plasma_raw

        # Radiative cooling
        eps_T = self.emissivity_T(emissivity_base, T_surface, composition)
        q_cool = self.q_surface_radiative_cooling(eps_T, T_surface, T_air)

        # Apply geometric factors to heating terms
        q_conv_effective = q_conv_total * aoa_factor * shape
        q_plasma_effective = q_plasma_raw * aoa_factor * shape
        q_fragment_effective = q_fragment * aoa_factor * shape
        q_catalytic_effective = q_catalytic * aoa_factor * shape

        # Total heating (before cooling)
        total_heating = q_conv_effective + q_plasma_effective + q_fragment_effective

        # Net heat flux (can be negative when cooling exceeds heating)
        q_net = total_heating - q_cool

        # Percentages of total heating (not net)
        if total_heating > 0:
            # Convective contribution (without catalytic, for breakdown)
            pct_conv = 100.0 * (q_conv_effective - q_catalytic_effective) / total_heating
            pct_rad = 100.0 * (q_plasma_effective + q_fragment_effective) / total_heating
            pct_cat = 100.0 * q_catalytic_effective / total_heating
        else:
            pct_conv = 0.0
            pct_rad = 0.0
            pct_cat = 0.0

        # Flow regime classification
        if blend_fm > 0.9:
            flow_regime = "free-molecular"
        elif blend_fm > 0.1:
            flow_regime = "transitional"
        else:
            flow_regime = "continuum"

        return HeatFluxBreakdown(
            q_convective=q_conv_effective - q_catalytic_effective,
            q_radiative_shock=q_plasma_effective + q_fragment_effective,
            q_radiative_cooling=q_cool,
            q_catalytic=q_catalytic_effective,
            q_fragment_reradiation=q_fragment_effective,
            q_net=q_net,
            pct_convective=pct_conv,
            pct_radiative=pct_rad,
            pct_catalytic=pct_cat,
            flow_regime=flow_regime,
            blend_factor=blend_fm,
        )

    def compute_biot_number(self, h_conv: float, radius: float, k_thermal: float) -> float:
        """
        Compute Biot number for thermal analysis.

        Bi = h·L / k

        where:
            h = convective heat transfer coefficient [W/(m²·K)]
            L = characteristic length (radius for sphere) [m]
            k = thermal conductivity [W/(m·K)]

        Bi < 0.1: Lumped capacitance valid (uniform temperature)
        Bi > 0.1: Temperature gradients are significant

        Parameters
        ----------
        h_conv : float
            Effective heat transfer coefficient [W/(m²·K)]
        radius : float
            Meteor radius [m]
        k_thermal : float
            Thermal conductivity [W/(m·K)]

        Returns
        -------
        float
            Biot number
        """
        return h_conv * radius / max(1e-9, k_thermal)

    def core_conduction_factor(self, radius: float, density: float, cp: float,
                                composition: str, dt: float) -> float:
        """
        Compute the fraction of surface heat that reaches the core per timestep.

        Uses thermal diffusivity to estimate heat penetration depth:
            δ ~ √(α·t) where α = k/(ρ·cp)

        The conduction factor increases for:
        - Higher thermal conductivity (metals)
        - Smaller radius (less distance to travel)
        - Longer timesteps

        Parameters
        ----------
        radius : float
            Meteor radius [m]
        density : float
            Material density [kg/m³]
        cp : float
            Specific heat [J/(kg·K)]
        composition : str
            Material type for thermal conductivity lookup
        dt : float
            Timestep [s]

        Returns
        -------
        float
            Conduction factor [0-1]
        """
        k_thermal = THERMAL_CONDUCTIVITY.get(composition, 2.0)
        alpha = k_thermal / max(1e-9, density * cp)  # Thermal diffusivity [m²/s]

        # Characteristic heat diffusion length scale
        diffusion_depth = math.sqrt(alpha * dt)

        # Fraction of radius that heat penetrates
        penetration_fraction = diffusion_depth / max(1e-6, radius)

        # For a sphere, the effective conduction factor accounts for
        # the decreasing area as heat moves inward
        # Use a saturating function to cap at reasonable values
        factor = min(0.5, penetration_fraction * 0.3)

        return factor

    def update_surface_and_core_temps(self, altitude_m: float, A: float, m: float, cp: float,
                                      emissivity: float, rho: float, v: float,
                                      T_surface: float, T_core: float, T_air: float,
                                      dt: float, composition: str, angle_rad: float,
                                      fragment_view_factor: float = 0.0,
                                      radius: float = None, density: float = None):
        """
        Update surface and core temperatures for one timestep.

        Uses a two-node thermal model:
        1. Surface node receives net heat flux and radiates
        2. Core node receives heat via conduction from surface

        The conduction rate depends on thermal diffusivity and body size.

        Parameters
        ----------
        altitude_m : float
            Altitude [m]
        A : float
            Cross-sectional area [m²]
        m : float
            Mass [kg]
        cp : float
            Specific heat [J/(kg·K)]
        emissivity : float
            Surface emissivity
        rho : float
            Atmospheric density [kg/m³]
        v : float
            Velocity [m/s]
        T_surface : float
            Current surface temperature [K]
        T_core : float
            Current core temperature [K]
        T_air : float
            Ambient air temperature [K]
        dt : float
            Timestep [s]
        composition : str
            Material type
        angle_rad : float
            Entry angle [rad]
        fragment_view_factor : float
            View factor for fragment re-radiation
        radius : float, optional
            Meteor radius [m] for conduction calculation
        density : float, optional
            Material density [kg/m³]

        Returns
        -------
        tuple
            (T_surface_new, T_core_new, P_net)
        """
        q_net = self.net_heat_flux(altitude_m, rho, v, emissivity, T_surface, T_air,
                                   composition, angle_rad, fragment_view_factor)
        P_net = q_net * max(1e-9, A)

        # Surface temperature change
        thermal_mass_surface = max(1e-9, m) * max(1e-9, cp)
        dTdt_surface = P_net / thermal_mass_surface

        # Conduction to core - use physics-based model if radius available
        if radius is not None and density is not None:
            conduction_factor = self.core_conduction_factor(radius, density, cp, composition, dt)
        else:
            # Fallback to simple fixed fraction
            conduction_factor = 0.10

        # Heat conducted to core from temperature difference
        k_thermal = THERMAL_CONDUCTIVITY.get(composition, 2.0)
        if radius is not None and radius > 1e-6:
            # Fourier's law approximation: Q = k·A·ΔT/L
            # For sphere, use effective area and path length
            delta_T = T_surface - T_core
            Q_conduction = k_thermal * A * delta_T / (0.5 * radius)  # W
            dTdt_core = Q_conduction / thermal_mass_surface
        else:
            # Simple fraction-based fallback
            dTdt_core = conduction_factor * dTdt_surface

        T_surface_new = T_surface + dTdt_surface * dt
        T_core_new = T_core + dTdt_core * dt

        # Ensure core doesn't exceed surface (thermodynamic consistency)
        if T_core_new > T_surface_new:
            T_core_new = T_surface_new

        return T_surface_new, T_core_new, P_net

    def ablation_rate(self, P_net: float, L_vap: float, T_surface: float,
                       fusion_T: float, vapor_T: float) -> float:
        """
        Compute mass ablation rate [kg/s].

        Ablation efficiency increases with temperature:
        - Below fusion: no ablation
        - Between fusion and vaporization: modest ablation (spallation, melting)
        - Above vaporization: vigorous ablation

        Parameters
        ----------
        P_net : float
            Net heating power [W]
        L_vap : float
            Latent heat of vaporization [J/kg]
        T_surface : float
            Surface temperature [K]
        fusion_T : float
            Fusion/melting temperature [K]
        vapor_T : float
            Vaporization temperature [K]

        Returns
        -------
        float
            Mass loss rate [kg/s]
        """
        if T_surface < 0.9 * fusion_T:
            return 0.0

        superheat = max(0.0, T_surface - fusion_T)
        if T_surface < vapor_T:
            # Below vaporization: melting and spallation
            f = 0.15 + 0.25 * (superheat / (superheat + 500.0))
        else:
            # Above vaporization: efficient mass loss
            f = 0.5 + 0.35 * (superheat / (superheat + 1000.0))

        f = min(0.9, f)
        return (f * max(0.0, P_net)) / max(1e-9, L_vap)
