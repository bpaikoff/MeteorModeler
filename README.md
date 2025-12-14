# MeteorModeler

A physics-based simulation of meteoroid atmospheric entry, ablation, and thermal evolution. This simulator models the descent of meteoroids through Earth's atmosphere, tracking temperature evolution, phase transitions, mass loss, and structural breakup.

## Overview

MeteorModeler implements a comprehensive physical model for meteor entry phenomena, including:

- **Aerodynamic Heating**: Convective and radiative heating using Sutton-Graves and shock-layer radiation models
- **Thermal Evolution**: Two-node thermal model with physics-based heat conduction
- **Phase Transitions**: Clausius-Clapeyron based vapor pressure with sublimation support
- **Ablation**: Mass loss through melting, vaporization, and spallation
- **Structural Breakup**: Weibull-scaled strength model with aerodynamic loading
- **Trajectory Dynamics**: Drag deceleration and gravitational acceleration
- **Plasma Physics**: Shock layer ionization, electron density, and plasma parameters
- **Real-time Visualization**: Heat flux breakdown and plasma stats with pause controls

## Physical Models

### Atmospheric Model

The atmosphere follows the U.S. Standard Atmosphere 1976 with extensions:

| Layer | Altitude (km) | Temperature Profile |
|-------|---------------|---------------------|
| Troposphere | 0-11 | 288.15 K, lapse rate -6.5 K/km |
| Stratosphere (lower) | 11-20 | Isothermal 216.65 K |
| Stratosphere (upper) | 20-32 | Lapse rate +1.0 K/km |
| Stratopause | 32-47 | Lapse rate +2.8 K/km |
| Mesosphere | 47-85 | Cooling to mesopause (~187 K) |
| Thermosphere | 85-700 | Exponential rise to ~1000 K |

Density follows an exponential profile with scale height 8.5 km.

### Heating Model

The net heat flux combines multiple mechanisms:

```
q_net = (q_convective + q_plasma) × f_angle × f_shape - q_radiative_cooling
```

**Convective Heating** (Sutton-Graves like):
```
q_conv ~ k × ρ^0.5 × v^3
```

**Shock-Layer Radiation**:
```
q_rad ~ k × ρ^0.8 × v^2.9
```

**Radiative Cooling** (Stefan-Boltzmann):
```
q_cool = ε × σ × (T_surface^4 - T_ambient^4)
```

The model transitions smoothly between free-molecular and continuum flow regimes based on altitude (Knudsen number proxy).

### Heat Flux Breakdown

The simulation tracks individual contributions to surface heating:

| Component | Description |
|-----------|-------------|
| Convective | Stagnation point convective heating |
| Radiative | Shock-layer thermal radiation |
| Catalytic | Surface catalytic recombination (~20% of convective) |
| Cooling | Stefan-Boltzmann surface radiation loss |
| Net | Total heat flux into the surface |

Contributions are displayed as percentages in a stacked bar visualization.

### Plasma Physics

The shock layer ahead of the meteor creates a high-temperature plasma. The model calculates:

**Shock Properties** (Rankine-Hugoniot relations):
```
T_shock = f(v, T_freestream)   # High-enthalpy correlation
ρ₂/ρ₁ = (γ+1)M² / (2 + (γ-1)M²)
δ/R ≈ 0.143 × exp(3.24/M²)    # Shock standoff
```

**Ionization** (Saha equation):
```
α²/(1-α) = K(T)/n    where K(T) ∝ T^(3/2) × exp(-χ/kT)
n_e = α × n_gas      # Electron density
```

**Plasma Parameters**:
- **Debye Length**: λ_D = √(ε₀kT_e / n_e·e²) — shielding distance
- **Plasma Frequency**: ω_pe = √(n_e·e² / ε₀·m_e) — oscillation frequency
- **Electrical Conductivity**: Spitzer formula with Coulomb logarithm
- **Magnetic Reynolds Number**: R_m = μ₀·σ·v·L — MHD coupling

**Dissociation**:
```
O₂ → 2O   (θ_d ≈ 59,500 K)
N₂ → 2N   (θ_d ≈ 113,000 K)
```

At typical meteor velocities (10-20 km/s), shock temperatures reach 5,000-15,000 K, causing significant dissociation and partial ionization.

### Thermal Conduction

Heat transfer from surface to core uses thermal diffusivity:

```
α = k / (ρ × cp)    [m²/s]
```

The penetration depth scales as:
```
δ ~ √(α × t)
```

This captures the physical reality that:
- Iron meteorites (high k) heat through quickly
- Stone/carbonaceous (low k) maintain cold cores longer
- Larger bodies take longer to heat through

### Phase Transitions

The phase model uses Clausius-Clapeyron vapor pressure:

```
P = P_ref × exp(-(L/R) × (1/T - 1/T_ref))
```

Key physics:
- **Above triple point pressure**: solid → liquid → vapor
- **Below triple point pressure**: solid → vapor (sublimation)
- **Energy-limited transitions**: Phase changes require both thermodynamic drive AND sufficient energy input

This is critical for ice meteoroids, which sublimate directly at high altitudes where ambient pressure is below 611 Pa. The model includes material-specific triple points and minimum sublimation temperatures to prevent instantaneous phase changes.

### Material Properties

| Property | Iron | Stone | Ice | Carbonaceous |
|----------|------|-------|-----|--------------|
| Density (kg/m³) | 7800 | 3500 | 917 | 2200 |
| Specific Heat (J/kg·K) | 450 | 800 | 2090 | 900 |
| Thermal Conductivity (W/m·K) | 80 | 2 | 2.2 | 1 |
| Melting Point (K) | 1809 | 1700 | 273 | 1500 |
| Vaporization Point (K) | 3134 | 2800 | 373 | 2200 |
| Tensile Strength (Pa) | 10⁸ | 2×10⁷ | 10⁶ | 10⁷ |

### Ablation

Mass loss rate depends on surface temperature relative to phase transition points:

```
ṁ = f(T) × P_net / L_vap
```

Where f(T) is an efficiency factor:
- Below 0.9 × T_fusion: no ablation
- Between fusion and vaporization: 15-40% efficiency (melting, spallation)
- Above vaporization: 50-85% efficiency

### Structural Breakup

The breakup model implements **progressive fragmentation** rather than binary failure:

1. **Weibull Size Scaling**: Larger bodies have lower effective strength
   ```
   σ_eff ~ σ_base × (r_ref / r)^(1/m)
   ```

2. **Phase Weakening**: Molten/vapor phases reduce strength significantly (10-60% of base)

3. **Thermal Stress**: Temperature gradients reduce structural integrity (25% per 1000K gradient)

4. **Aerodynamic Loading**: Dynamic pressure and bending moments

**Progressive Fragmentation**:
- **Onset threshold** (stress ratio > 0.7): Fragment shedding begins
- **Transition zone** (0.7 - 1.5): Probabilistic fragmentation, rate increases with stress
- **Cascade threshold** (stress ratio > 1.5): Rapid breakup with 50-90% mass loss

Fragment masses follow a power-law distribution (α = 1.8), consistent with observed meteoroid fragmentation.

## Installation

### Requirements

- Python 3.10+
- pygame (for visualization)
- matplotlib (optional, for analysis plots)

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/MeteorModeler.git
cd MeteorModeler

# Install dependencies
pip install pygame matplotlib
```

## Usage

### Interactive Simulation

```bash
python main.py
```

The simulator will prompt you to select:
1. **Composition**: iron, stone, ice, or carbonaceous
2. **Diameter**: initial meteor diameter in meters
3. **Lifeform**: hypothetical organism type for survival analysis

### Keyboard Controls

| Key | Action |
|-----|--------|
| **Space** or **P** | Pause/Resume simulation |
| **Escape** | Quit simulation |

### User Interface

The simulation window displays three panels:

- **Left Panel**: Meteor status (altitude, velocity, temperature, mass, phase)
- **Center**: Visualization of meteor descent with thermal trail
- **Right Panel**: Real-time physics data
  - *Heat Flux Breakdown*: Individual heating contributions with percentage bar
  - *Plasma Properties*: Shock layer, ionization, and flow parameters

### Survival Analysis

The survival model tracks thermal exposure throughout the meteor's descent:
- If the meteor **fully ablates** before reaching the ground, survival is **0%**
- If the meteor **reaches the ground**, survival probability is computed using an exponential hazard model based on cumulative thermal exposure above the organism's critical temperature

This reflects the physical reality that biological material cannot survive if its carrier is completely destroyed.

### Validation/Headless Mode

```bash
python validate.py
```

Runs a simulation and outputs trajectory data without visualization.

### Monte Carlo Analysis

```bash
python montecarlo.py
```

Runs multiple simulations with varied parameters to analyze statistical outcomes.

## Project Structure

```
MeteorModeler/
├── main.py              # Main entry point with UI
├── equations.py         # Coordinator for physics models
├── lifeforms.py         # Survival probability model
├── ui.py                # Pygame visualization
├── validate.py          # Headless validation runs
├── montecarlo.py        # Statistical analysis
│
├── physics/             # Physics simulation modules
│   ├── __init__.py      # Package exports
│   ├── atmosphere.py    # Atmospheric density, temperature, pressure
│   ├── trajectory.py    # Trajectory dynamics and drag
│   ├── heating.py       # Aerodynamic and radiative heating
│   ├── plasma.py        # Shock layer plasma physics
│   ├── ablation.py      # Mass loss calculations
│   ├── phase.py         # Phase transition logic
│   ├── vapor_pressure.py    # Clausius-Clapeyron vapor pressure
│   ├── energy_balance.py    # Energy partitioning model
│   ├── breakup.py       # Structural failure and fragmentation
│   └── vaporization_detector.py  # Sustained vaporization detection
│
└── entities/            # Data structures and simulation objects
    ├── __init__.py      # Package exports
    ├── meteor.py        # Meteor object with state tracking
    ├── fragment.py      # Fragment thermal evolution
    ├── particles.py     # Particle/fragment data class
    └── compositions.py  # Material property database
```

## Simulation Parameters

### Entry Conditions (defaults)

| Parameter | Value | Notes |
|-----------|-------|-------|
| Entry altitude | 120 km | Top of thermosphere |
| Entry velocity | 19 km/s | Typical NEA velocity |
| Entry angle | 45° | From vertical |
| Timestep | 0.02 s | For numerical stability |

### Physical Constants

| Constant | Value | Symbol |
|----------|-------|--------|
| Stefan-Boltzmann | 5.67×10⁻⁸ W/(m²·K⁴) | σ |
| Universal gas constant | 8.314 J/(mol·K) | R |
| Sea-level density | 1.225 kg/m³ | ρ₀ |
| Scale height | 8500 m | H |
| Sea-level pressure | 101325 Pa | P₀ |

## Scientific References

1. **Meteor Physics**: Ceplecha, Z. et al. (1998). "Meteor Phenomena and Bodies." Space Science Reviews 84, 327-471.

2. **Atmospheric Entry Heating**: Sutton, K. & Graves, R.A. (1971). "A General Stagnation-Point Convective-Heating Equation for Arbitrary Gas Mixtures." NASA TR R-376.

3. **Vapor Pressure**: Wagner, W. & Pruss, A. (2002). "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance." J. Phys. Chem. Ref. Data 31, 387-535.

4. **Meteoroid Strength**: Popova, O.P. et al. (2011). "Very Low Strengths of Interplanetary Meteoroids and Small Asteroids." Meteoritics & Planetary Science 46, 1525-1550.

5. **U.S. Standard Atmosphere**: NOAA/NASA/USAF (1976). U.S. Standard Atmosphere, 1976.

6. **Plasma Physics**: Park, C. (1989). "Nonequilibrium Hypersonic Aerothermodynamics." Wiley-Interscience.

7. **High-Temperature Gas Dynamics**: Anderson, J.D. (2006). "Hypersonic and High-Temperature Gas Dynamics." AIAA Education Series.

## Limitations

- **1D Thermal Model**: Uses a two-node (surface/core) approximation rather than full 3D heat equation
- **Simplified Aerodynamics**: No full CFD; uses empirical correlations
- **Spherical Geometry**: Assumes spherical meteoroid shape
- **Single-Component Materials**: Does not model heterogeneous composition
- **No Atmospheric Winds**: Assumes still atmosphere
- **Simplified Fragmentation**: Power-law mass distribution, no detailed fracture mechanics

## Contributing

Contributions are welcome! Areas for improvement include:

- 3D thermal conduction (finite element)
- Non-spherical shape effects
- Atmospheric wind/turbulence
- Multi-component ablation (differential volatility)
- Real meteor light curve comparison
- GPU acceleration for Monte Carlo

## License

MIT License - see LICENSE file for details.

## Acknowledgments

This simulator was developed for educational purposes to demonstrate the physics of meteor entry. The models are simplified representations suitable for understanding fundamental principles but should not be used for spacecraft reentry design or impact hazard assessment.
