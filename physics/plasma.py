# plasma.py
"""
Plasma physics model for meteor entry shock layer.

This module calculates plasma properties in the high-temperature shock layer
formed during hypersonic atmospheric entry:

1. Shock temperature via Rankine-Hugoniot relations
2. Ionization fraction via Saha equation
3. Electron/ion number densities
4. Electrical conductivity (Spitzer formula)
5. Debye length and plasma frequency
6. Magnetic Reynolds number for MHD effects
7. Radiative properties

References:
- Park, C. (1989) "Nonequilibrium Hypersonic Aerothermodynamics"
- Panesi, M. et al. (2009) "Fire II flight experiment analysis"
- Anderson, J.D. (2006) "Hypersonic and High-Temperature Gas Dynamics"
"""
import math
from dataclasses import dataclass


# Physical constants
K_BOLTZMANN = 1.380649e-23      # J/K
E_CHARGE = 1.602176634e-19      # C
EPSILON_0 = 8.854187817e-12     # F/m (vacuum permittivity)
M_ELECTRON = 9.10938e-31        # kg
M_PROTON = 1.67262e-27          # kg
AVOGADRO = 6.02214076e23        # mol^-1
R_GAS = 8.314462                # J/(mol·K)
PLANCK = 6.62607e-34            # J·s
C_LIGHT = 2.998e8               # m/s

# Air composition (simplified)
MEAN_MOLECULAR_MASS_AIR = 28.97e-3 / AVOGADRO  # kg per molecule (~4.81e-26 kg)
GAMMA_AIR_COLD = 1.4            # Ratio of specific heats (cold air)
GAMMA_AIR_HOT = 1.2             # At high T, gamma decreases due to vibrational modes


@dataclass
class PlasmaStats:
    """Container for plasma properties in the shock layer."""
    # Shock properties
    shock_standoff: float = 0.0         # Shock standoff distance [m]
    shock_temperature: float = 0.0       # Post-shock temperature [K]
    shock_density_ratio: float = 1.0     # ρ2/ρ1 across shock
    mach_number: float = 0.0             # Freestream Mach number

    # Ionization
    ionization_fraction: float = 0.0     # α = n_e / n_total
    electron_density: float = 0.0        # n_e [m^-3]
    ion_density: float = 0.0             # n_i [m^-3]
    mean_charge_state: float = 0.0       # Average Z for ions

    # Plasma parameters
    debye_length: float = 0.0            # λ_D [m]
    plasma_frequency: float = 0.0        # ω_pe [rad/s]
    plasma_parameter: float = 0.0        # Λ = n_e * λ_D^3 (coupling)
    electron_temperature: float = 0.0    # T_e [K] (can differ from T_heavy)

    # Transport
    electrical_conductivity: float = 0.0  # σ [S/m]
    thermal_conductivity_plasma: float = 0.0  # k_plasma [W/(m·K)]
    magnetic_reynolds: float = 0.0        # Rm for MHD effects

    # Radiative
    radiative_power_density: float = 0.0  # P_rad [W/m³]
    optical_depth: float = 0.0            # τ (dimensionless)
    planck_mean_absorption: float = 0.0   # κ_P [m^-1]

    # Species
    dissociation_N2: float = 0.0         # Fraction of N2 dissociated
    dissociation_O2: float = 0.0         # Fraction of O2 dissociated

    # Derived quantities
    reynolds_number: float = 0.0         # Re = ρvL/μ
    knudsen_number: float = 0.0          # Kn = λ/L (mean free path / size)


@dataclass
class HeatFluxBreakdown:
    """Breakdown of individual heat flux contributions."""
    q_convective: float = 0.0            # Convective heating [W/m²]
    q_radiative_shock: float = 0.0       # Shock layer radiation [W/m²]
    q_radiative_cooling: float = 0.0     # Surface radiative cooling [W/m²]
    q_catalytic: float = 0.0             # Catalytic recombination heating [W/m²]
    q_fragment_reradiation: float = 0.0  # Re-radiation from fragments [W/m²]
    q_net: float = 0.0                   # Net heat flux [W/m²]

    # Percentages (of total heating, not net)
    pct_convective: float = 0.0
    pct_radiative: float = 0.0
    pct_catalytic: float = 0.0

    # Flow regime
    flow_regime: str = "continuum"       # "free-molecular", "transitional", "continuum"
    blend_factor: float = 0.0            # 0 = continuum, 1 = free-molecular


class PlasmaModel:
    """
    Calculates plasma properties in the meteor shock layer.

    The shock layer forms when the meteor travels supersonically through
    the atmosphere. At hypersonic speeds (M > 5), the shocked gas reaches
    temperatures high enough to cause:
    1. Molecular dissociation (N2, O2 → N, O)
    2. Electronic excitation
    3. Ionization (N, O → N+, O+ + e-)

    These effects create a partially ionized plasma that radiates intensely.
    """

    def __init__(self):
        # Ionization energies [eV]
        self.ionization_energy = {
            "N": 14.534,   # Nitrogen
            "O": 13.618,   # Oxygen
            "N2": 15.58,   # Molecular nitrogen (first ionization)
            "O2": 12.07,   # Molecular oxygen
            "air": 14.0,   # Effective average for air
        }

        # Dissociation energies [eV]
        self.dissociation_energy = {
            "N2": 9.79,    # N2 → 2N
            "O2": 5.12,    # O2 → 2O
        }

        # Characteristic dissociation temperatures [K]
        self.dissociation_temp = {
            "N2": 113_000,  # θ_d for N2
            "O2": 59_500,   # θ_d for O2
        }

    def compute_mach_number(self, velocity: float, temperature: float,
                            gamma: float = GAMMA_AIR_COLD) -> float:
        """
        Compute Mach number from velocity and atmospheric temperature.

        M = v / a, where a = sqrt(γRT/M_mol)
        """
        # Speed of sound in air
        R_specific = R_GAS / 0.02897  # J/(kg·K) for air
        a = math.sqrt(gamma * R_specific * max(100.0, temperature))
        return velocity / max(1.0, a)

    def shock_temperature_rankine_hugoniot(self, v: float, T_freestream: float,
                                           gamma: float = GAMMA_AIR_COLD) -> float:
        """
        Post-shock temperature from normal shock relations.

        For a strong shock (M >> 1):
            T2/T1 ≈ (2γ(γ-1)M²) / (γ+1)²

        This uses the calorically perfect gas assumption which breaks down
        at very high temperatures due to real gas effects.
        """
        M = self.compute_mach_number(v, T_freestream, gamma)

        if M < 1.0:
            return T_freestream

        M2 = M * M
        gp1 = gamma + 1.0
        gm1 = gamma - 1.0

        # Normal shock temperature ratio
        numerator = (1.0 + 0.5 * gm1 * M2) * (2.0 * gamma * M2 / gm1 - 1.0)
        denominator = M2 * (2.0 * gamma / gm1 + 0.5 * gm1)

        T_ratio = numerator / max(1e-9, denominator)
        T_shock = T_freestream * max(1.0, T_ratio)

        # Cap at physically reasonable values (~50,000 K for meteor entry)
        return min(50_000.0, T_shock)

    def shock_temperature_high_enthalpy(self, v: float, T_freestream: float) -> float:
        """
        Post-shock temperature accounting for real gas effects.

        At high temperatures (> 2000 K), vibrational modes become excited
        and molecules begin to dissociate, which absorbs energy and limits
        temperature rise. This uses an empirical correlation.

        For meteor velocities (10-70 km/s), equilibrium shock temperatures
        range from ~5,000 K to ~15,000 K due to these effects.
        """
        # Stagnation enthalpy: h0 = v²/2 + cp*T
        # For high-enthalpy flows, use effective cp that accounts for
        # dissociation and ionization energy sinks

        v_km_s = v / 1000.0

        # Empirical correlation for equilibrium shock temperature
        # Based on Park (1989) and flight data
        if v_km_s < 3.0:
            # Low velocity: nearly ideal gas behavior
            return self.shock_temperature_rankine_hugoniot(v, T_freestream)
        elif v_km_s < 8.0:
            # Moderate hypersonic: dissociation effects
            T_shock = 2500.0 + 900.0 * v_km_s
        elif v_km_s < 15.0:
            # High hypersonic: strong dissociation and some ionization
            T_shock = 6000.0 + 400.0 * v_km_s
        else:
            # Very high velocity (meteor regime): ionization plateau
            T_shock = 10000.0 + 150.0 * v_km_s

        return min(50_000.0, T_shock)

    def shock_density_ratio(self, mach: float, gamma: float = GAMMA_AIR_COLD) -> float:
        """
        Density ratio across a normal shock.

        ρ2/ρ1 = (γ+1)M² / (2 + (γ-1)M²)

        For strong shocks, this approaches (γ+1)/(γ-1) ≈ 6 for air.
        """
        if mach < 1.0:
            return 1.0

        M2 = mach * mach
        return (gamma + 1.0) * M2 / (2.0 + (gamma - 1.0) * M2)

    def shock_standoff_distance(self, radius: float, mach: float,
                                 density_ratio: float) -> float:
        """
        Shock standoff distance for a blunt body.

        δ/R ≈ 0.143 * exp(3.24/M²) for high Mach

        A simpler correlation: δ/R ≈ ρ1/ρ2 (inverse density ratio)
        """
        if mach < 1.5:
            return 0.5 * radius  # Weak shock

        # Empirical relation for blunt bodies
        delta_over_R = 0.143 * math.exp(3.24 / (mach * mach))

        # Also scale with density ratio
        delta_over_R = min(delta_over_R, 1.0 / max(1.0, density_ratio))

        return delta_over_R * radius

    def dissociation_fraction(self, temperature: float, species: str = "O2") -> float:
        """
        Equilibrium dissociation fraction using simplified model.

        Uses characteristic dissociation temperature:
            α_d ≈ 1 - exp(-T / θ_d)

        More accurate: Solve coupled equilibrium equations
        """
        theta_d = self.dissociation_temp.get(species, 60_000)

        if temperature < 2000.0:
            return 0.0

        # Simplified equilibrium dissociation
        alpha = 1.0 - math.exp(-temperature / theta_d)
        return min(1.0, max(0.0, alpha))

    def saha_ionization_fraction(self, temperature: float,
                                  number_density: float,
                                  species: str = "air") -> float:
        """
        Ionization fraction from the Saha equation.

        For equilibrium ionization:
            n_e * n_i / n_0 = (2/Λ³) * (g_i/g_0) * exp(-χ/kT)

        where Λ = h/√(2πm_e kT) is the thermal de Broglie wavelength.

        Simplifying for singly-ionized species with charge neutrality (n_e = n_i):
            α² / (1-α) = K(T) / n_total
        """
        chi_eV = self.ionization_energy.get(species, 14.0)
        chi_J = chi_eV * E_CHARGE  # Convert to Joules

        if temperature < 3000.0:
            return 0.0

        # Thermal de Broglie wavelength
        lambda_db = PLANCK / math.sqrt(2.0 * math.pi * M_ELECTRON * K_BOLTZMANN * temperature)

        # Saha constant (simplified, assuming g_i/g_0 ≈ 1)
        K_saha = (2.0 / (lambda_db ** 3)) * math.exp(-chi_J / (K_BOLTZMANN * temperature))

        # Solve quadratic: α²·n + α·K - K = 0 where n = number density
        # Actually: n_e² / (n - n_e) = K, so n_e² + K·n_e - K·n = 0
        # Let x = n_e/n = α: x²·n + K·x - K = 0...
        # Simpler: x² / (1-x) = K/n

        ratio = K_saha / max(1e10, number_density)

        # Quadratic solution for α: α² + ratio·α - ratio = 0
        # α = (-ratio + sqrt(ratio² + 4·ratio)) / 2
        discriminant = ratio * ratio + 4.0 * ratio
        alpha = (-ratio + math.sqrt(discriminant)) / 2.0

        return min(1.0, max(0.0, alpha))

    def electron_density(self, ionization_fraction: float,
                          gas_number_density: float) -> float:
        """
        Electron number density [m^-3].

        n_e = α * n_gas
        """
        return ionization_fraction * gas_number_density

    def debye_length(self, electron_density: float,
                      electron_temperature: float) -> float:
        """
        Debye length - characteristic shielding distance in plasma.

        λ_D = √(ε₀ k T_e / (n_e e²))

        Typically ~10⁻⁶ to 10⁻⁴ m in meteor plasmas.
        """
        if electron_density < 1e10 or electron_temperature < 1000.0:
            return float('inf')

        numerator = EPSILON_0 * K_BOLTZMANN * electron_temperature
        denominator = electron_density * E_CHARGE * E_CHARGE

        return math.sqrt(numerator / denominator)

    def plasma_frequency(self, electron_density: float) -> float:
        """
        Electron plasma frequency [rad/s].

        ω_pe = √(n_e e² / (ε₀ m_e))

        This is the characteristic frequency of plasma oscillations.
        For meteor plasmas: ~10⁹ to 10¹¹ rad/s
        """
        if electron_density < 1e10:
            return 0.0

        return math.sqrt(electron_density * E_CHARGE**2 / (EPSILON_0 * M_ELECTRON))

    def plasma_parameter(self, electron_density: float, debye_length: float) -> float:
        """
        Plasma coupling parameter Λ = n_e · λ_D³.

        Λ >> 1: Weakly coupled (ideal) plasma
        Λ ~ 1:  Strongly coupled plasma

        For meteor entry plasmas, typically Λ ~ 10³-10⁶ (weakly coupled).
        """
        if debye_length == float('inf') or electron_density < 1e10:
            return float('inf')

        return electron_density * (debye_length ** 3)

    def electrical_conductivity_spitzer(self, electron_temperature: float,
                                         electron_density: float,
                                         ion_charge: float = 1.0) -> float:
        """
        Electrical conductivity from Spitzer formula.

        σ = (3/√(2π)) · (4πε₀)² · (kT_e)^(3/2) / (Z e² m_e^(1/2) · ln(Λ))

        For partially ionized gases, need to also consider electron-neutral
        collisions which can dominate at low ionization fractions.

        Returns conductivity in S/m.
        """
        if electron_temperature < 3000.0 or electron_density < 1e10:
            return 0.0

        # Coulomb logarithm (typical range 5-20 for plasmas)
        lambda_D = self.debye_length(electron_density, electron_temperature)
        if lambda_D == float('inf'):
            return 0.0

        # Approximate mean interparticle distance
        b_min = E_CHARGE**2 / (4.0 * math.pi * EPSILON_0 * K_BOLTZMANN * electron_temperature)

        if lambda_D > b_min:
            ln_coulomb = math.log(lambda_D / b_min)
        else:
            ln_coulomb = 5.0  # Minimum reasonable value

        ln_coulomb = max(5.0, min(25.0, ln_coulomb))

        # Spitzer conductivity
        factor = 3.0 / math.sqrt(2.0 * math.pi)
        eps_factor = (4.0 * math.pi * EPSILON_0) ** 2
        kT_32 = (K_BOLTZMANN * electron_temperature) ** 1.5
        denominator = ion_charge * E_CHARGE**2 * math.sqrt(M_ELECTRON) * ln_coulomb

        sigma = factor * eps_factor * kT_32 / denominator

        # For partially ionized plasma, reduce by ionization fraction effect
        # (electron-neutral collisions become important)
        return sigma

    def thermal_conductivity_plasma(self, electron_temperature: float,
                                     electron_density: float) -> float:
        """
        Thermal conductivity of the plasma [W/(m·K)].

        For fully ionized plasma, thermal conductivity is dominated by
        electrons and scales as T^(5/2).
        """
        sigma = self.electrical_conductivity_spitzer(electron_temperature, electron_density)

        if sigma < 1e-6:
            return 0.0

        # Wiedemann-Franz law relation
        # κ = L · σ · T, where L ≈ 2.44e-8 W·Ω/K² (Lorenz number)
        L_lorenz = 2.44e-8

        return L_lorenz * sigma * electron_temperature

    def magnetic_reynolds_number(self, velocity: float, length_scale: float,
                                  conductivity: float) -> float:
        """
        Magnetic Reynolds number for MHD effects.

        Rm = μ₀ · σ · v · L

        Rm >> 1: Magnetic field is frozen into the plasma
        Rm << 1: Magnetic field diffuses freely

        For meteor entry, Rm is typically small (~0.01-1).
        """
        MU_0 = 4.0 * math.pi * 1e-7  # H/m
        return MU_0 * conductivity * velocity * length_scale

    def radiative_power_density(self, temperature: float,
                                 number_density: float,
                                 ionization_fraction: float) -> float:
        """
        Volumetric radiative power from the shock layer [W/m³].

        This includes:
        - Bound-bound transitions (line radiation)
        - Bound-free transitions (recombination radiation)
        - Free-free transitions (bremsstrahlung)

        Simplified model based on Park's correlations.
        """
        if temperature < 3000.0:
            return 0.0

        # Bremsstrahlung (free-free) power
        # P_ff ∝ n_e² · T^(1/2)
        n_e = ionization_fraction * number_density

        if n_e < 1e10:
            return 0.0

        # Gaunt factor (quantum correction) ~ 1.2 for average
        g_ff = 1.2

        # Bremsstrahlung emission coefficient
        # P_ff = 1.42e-40 · g_ff · n_e · n_i · T^(1/2) [W/m³] (CGS-based formula converted)
        P_ff = 1.42e-40 * g_ff * n_e * n_e * math.sqrt(temperature)

        # Add line radiation contribution (can be 10-100x bremsstrahlung)
        # Simplified: use temperature-dependent multiplier
        if temperature < 6000.0:
            line_multiplier = 5.0
        elif temperature < 10000.0:
            line_multiplier = 20.0
        elif temperature < 15000.0:
            line_multiplier = 50.0  # Strong N, O lines
        else:
            line_multiplier = 30.0  # Some lines saturate

        return P_ff * line_multiplier

    def optical_depth(self, absorption_coeff: float, path_length: float) -> float:
        """
        Optical depth τ = κ · L.

        τ << 1: Optically thin (radiation escapes freely)
        τ >> 1: Optically thick (radiation trapped)
        """
        return absorption_coeff * path_length

    def planck_mean_absorption(self, temperature: float,
                                number_density: float,
                                ionization_fraction: float) -> float:
        """
        Planck-mean absorption coefficient [m^-1].

        Approximate model for air plasma opacity.
        """
        if temperature < 4000.0:
            return 0.0

        n_e = ionization_fraction * number_density

        if n_e < 1e10:
            return 0.0

        # Simplified opacity model
        # κ ∝ n_e² / T^(7/2) for free-free
        kappa_ff = 1e-50 * n_e * n_e / (temperature ** 3.5)

        # Add bound-free contribution
        kappa_bf = 2.0 * kappa_ff  # Rough approximation

        return kappa_ff + kappa_bf

    def knudsen_number(self, mean_free_path: float, characteristic_length: float) -> float:
        """
        Knudsen number Kn = λ / L.

        Kn < 0.01: Continuum flow
        0.01 < Kn < 0.1: Slip flow
        0.1 < Kn < 10: Transitional
        Kn > 10: Free molecular flow
        """
        return mean_free_path / max(1e-9, characteristic_length)

    def mean_free_path(self, number_density: float, temperature: float) -> float:
        """
        Mean free path for gas molecules [m].

        λ = 1 / (√2 · n · σ)

        where σ ≈ πd² is the collision cross-section (d ~ 3.7e-10 m for air).
        """
        if number_density < 1e10:
            return float('inf')

        d_mol = 3.7e-10  # m, effective molecular diameter for air
        sigma_collision = math.pi * d_mol * d_mol

        return 1.0 / (math.sqrt(2.0) * number_density * sigma_collision)

    def reynolds_number(self, density: float, velocity: float,
                         length: float, temperature: float) -> float:
        """
        Reynolds number Re = ρvL/μ.

        Uses Sutherland's formula for air viscosity.
        """
        # Sutherland's formula for air
        mu_ref = 1.716e-5  # Pa·s at T_ref
        T_ref = 273.15     # K
        S = 110.4          # K (Sutherland constant)

        # At high temperatures, use power law
        if temperature < 500:
            mu = mu_ref * (temperature / T_ref)**1.5 * (T_ref + S) / (temperature + S)
        else:
            # High temperature: μ ∝ T^0.7 approximately
            mu = mu_ref * (temperature / T_ref)**0.7

        return density * velocity * length / max(1e-12, mu)

    def compute_all(self, velocity: float, altitude_m: float,
                     density: float, temperature_freestream: float,
                     meteor_radius: float) -> PlasmaStats:
        """
        Compute all plasma properties for current conditions.

        Parameters
        ----------
        velocity : float
            Meteor velocity [m/s]
        altitude_m : float
            Altitude [m]
        density : float
            Atmospheric density [kg/m³]
        temperature_freestream : float
            Freestream (ambient) temperature [K]
        meteor_radius : float
            Meteor radius [m]

        Returns
        -------
        PlasmaStats
            Complete plasma state
        """
        stats = PlasmaStats()

        # Mach number
        stats.mach_number = self.compute_mach_number(velocity, temperature_freestream)

        # Shock temperature (use high-enthalpy model for accuracy)
        stats.shock_temperature = self.shock_temperature_high_enthalpy(
            velocity, temperature_freestream
        )

        # Density ratio and shock standoff
        stats.shock_density_ratio = self.shock_density_ratio(stats.mach_number)
        stats.shock_standoff = self.shock_standoff_distance(
            meteor_radius, stats.mach_number, stats.shock_density_ratio
        )

        # Post-shock number density
        n_freestream = density / MEAN_MOLECULAR_MASS_AIR
        n_shock = n_freestream * stats.shock_density_ratio

        # Dissociation
        stats.dissociation_O2 = self.dissociation_fraction(
            stats.shock_temperature, "O2"
        )
        stats.dissociation_N2 = self.dissociation_fraction(
            stats.shock_temperature, "N2"
        )

        # Ionization (Saha equation)
        stats.ionization_fraction = self.saha_ionization_fraction(
            stats.shock_temperature, n_shock, "air"
        )

        # Electron density and temperature
        stats.electron_density = self.electron_density(
            stats.ionization_fraction, n_shock
        )
        stats.ion_density = stats.electron_density  # Charge neutrality
        stats.electron_temperature = stats.shock_temperature  # Thermal equilibrium assumption
        stats.mean_charge_state = 1.0  # Singly ionized approximation

        # Plasma parameters
        stats.debye_length = self.debye_length(
            stats.electron_density, stats.electron_temperature
        )
        stats.plasma_frequency = self.plasma_frequency(stats.electron_density)
        stats.plasma_parameter = self.plasma_parameter(
            stats.electron_density, stats.debye_length
        )

        # Transport properties
        stats.electrical_conductivity = self.electrical_conductivity_spitzer(
            stats.electron_temperature, stats.electron_density
        )
        stats.thermal_conductivity_plasma = self.thermal_conductivity_plasma(
            stats.electron_temperature, stats.electron_density
        )
        stats.magnetic_reynolds = self.magnetic_reynolds_number(
            velocity, meteor_radius, stats.electrical_conductivity
        )

        # Radiative properties
        stats.radiative_power_density = self.radiative_power_density(
            stats.shock_temperature, n_shock, stats.ionization_fraction
        )
        stats.planck_mean_absorption = self.planck_mean_absorption(
            stats.shock_temperature, n_shock, stats.ionization_fraction
        )
        stats.optical_depth = self.optical_depth(
            stats.planck_mean_absorption, stats.shock_standoff
        )

        # Flow characterization
        mfp = self.mean_free_path(n_freestream, temperature_freestream)
        stats.knudsen_number = self.knudsen_number(mfp, 2.0 * meteor_radius)
        stats.reynolds_number = self.reynolds_number(
            density, velocity, 2.0 * meteor_radius, temperature_freestream
        )

        return stats
