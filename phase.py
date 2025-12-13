# phase.py
"""
Phase transition model for meteoroid materials.

This module determines the thermodynamic phase of the material based on:
- Temperature relative to fusion and vaporization points
- Ambient pressure relative to vapor pressure (Clausius-Clapeyron)
- Net heating rate (energy-limited phase transitions)

Key physics:
- At low pressures (upper atmosphere), materials can sublimate directly
  from solid to vapor, bypassing the liquid phase
- BUT sublimation rate is energy-limited, not instantaneous
- The triple point defines where solid, liquid, and vapor coexist
- Phase transitions require both thermodynamic drive AND energy supply
"""


class PhaseModel:
    def __init__(self, vapor_model, detector):
        self.vapor = vapor_model
        self.detector = detector

        # Triple point temperatures [K] and pressures [Pa]
        self.triple_point = {
            "ice": {"T": 273.16, "P": 611.657},
            "iron": {"T": 1809.0, "P": 1.0},      # Iron: no practical sublimation
            "stone": {"T": 1700.0, "P": 1.0},     # Silicates: congruent melting
            "carbonaceous": {"T": 1500.0, "P": 10.0}
        }

        # Minimum temperature for significant sublimation rate [K]
        # Below this, sublimation is negligible even if thermodynamically favored
        self.min_sublimation_T = {
            "ice": 150.0,          # Ice starts sublimating noticeably above ~150K
            "iron": 2000.0,        # Iron needs extreme temperatures
            "stone": 1500.0,       # Silicates need high T
            "carbonaceous": 400.0  # Volatiles come off at lower T
        }

    def evaluate(self, mass, props, T_surface, ambient_pressure, material_key: str,
                 latent_fusion_remaining: float, latent_vapor_remaining: float,
                 P_net: float = 0.0):
        """
        Determine the current thermodynamic phase.

        The phase indicates what physical process dominates, but actual mass loss
        is controlled by energy balance, not phase alone.

        Parameters
        ----------
        mass : float
            Current mass [kg]
        props : dict
            Material properties including fusion_T, vapor_T, L_fusion
        T_surface : float
            Surface temperature [K]
        ambient_pressure : float
            Ambient atmospheric pressure [Pa]
        material_key : str
            Material type identifier
        latent_fusion_remaining : float
            Remaining fusion latent heat inventory [J]
        latent_vapor_remaining : float
            Remaining vaporization latent heat inventory [J]
        P_net : float
            Net heating power [W] - sublimation only active if P_net > 0

        Returns
        -------
        str
            Phase: "solid", "fusion", "liquid", "sublimation", or "vapor"
        """
        fusion_T = props["fusion_T"]
        vapor_T = props.get("vapor_T", fusion_T + 500.0)

        # Get material-specific thresholds
        tp = self.triple_point.get(material_key, {"T": fusion_T, "P": 100.0})
        T_triple = tp["T"]
        P_triple = tp["P"]
        T_min_sub = self.min_sublimation_T.get(material_key, 0.5 * fusion_T)

        # Calculate vapor pressure at current temperature
        Pvap = self.vapor.vapor_pressure(material_key, T_surface)

        # Calculate fusion fraction consumed (with safe division)
        L_fusion_total = mass * props["L_fusion"]
        if L_fusion_total > 1e-9:
            fusion_fraction_consumed = 1.0 - (latent_fusion_remaining / L_fusion_total)
        else:
            fusion_fraction_consumed = 1.0

        # Determine if sublimation/vaporization is thermodynamically active
        # Requires: 1) Pvap > Pamb, 2) T above minimum threshold, 3) Energy available
        sublimation_active = (
            Pvap > ambient_pressure and
            T_surface >= T_min_sub and
            P_net > 0  # Must have net energy input
        )

        # Phase determination logic
        if T_surface < fusion_T:
            # Below melting point
            if sublimation_active and ambient_pressure < P_triple:
                phase = "sublimation"
            else:
                phase = "solid"

        elif T_surface < vapor_T:
            # Between fusion and vaporization temperatures
            if ambient_pressure < P_triple and sublimation_active:
                # Below triple point with energy: direct sublimation
                phase = "sublimation"
            elif latent_fusion_remaining > 1e-6 and fusion_fraction_consumed < 0.9:
                # Still consuming fusion latent heat
                phase = "fusion"
            elif sublimation_active:
                # Above fusion, energy available: evaporating
                phase = "vapor"
            else:
                # Liquid phase stable
                phase = "liquid"

        else:
            # Above vaporization temperature
            if Pvap > ambient_pressure and P_net > 0:
                phase = "vapor"
            elif Pvap > ambient_pressure:
                # Thermodynamically wants to vaporize but no energy
                phase = "liquid"  # Treat as liquid, rate-limited
            else:
                # High pressure suppresses vaporization
                phase = "liquid"

        self.detector.update(phase, T_surface, fusion_T)
        return phase

    def is_sublimating(self, phase: str) -> bool:
        """Check if material is undergoing sublimation."""
        return phase in ("sublimation", "vapor")

    def core_sublimated(self) -> bool:
        """
        Check if the core has been destroyed by sustained high-temperature vaporization.

        This requires both sustained vapor phase AND high temperatures.
        """
        return self.detector.sustained_destruction()
