# phase.py
"""
Phase transition model for meteoroid materials.

This module determines the thermodynamic phase of the material based on:
- Temperature relative to fusion and vaporization points
- Ambient pressure relative to vapor pressure (Clausius-Clapeyron)
- Remaining latent heat inventories

Key physics:
- At low pressures (upper atmosphere), materials can sublimate directly
  from solid to vapor, bypassing the liquid phase
- The triple point defines where solid, liquid, and vapor coexist
- Above the triple point pressure: solid → liquid → vapor
- Below the triple point pressure: solid → vapor (sublimation)
"""


class PhaseModel:
    def __init__(self, vapor_model, detector):
        self.vapor = vapor_model
        self.detector = detector

        # Triple point pressures for phase path determination [Pa]
        self.triple_point_pressure = {
            "ice": 611.657,       # Water triple point
            "iron": 1.0,          # Iron sublimes only in extreme vacuum
            "stone": 1.0,         # Silicates effectively don't have liquid phase in vacuum
            "carbonaceous": 10.0  # Low due to volatile organics
        }

    def evaluate(self, mass, props, T_surface, ambient_pressure, material_key: str,
                 latent_fusion_remaining: float, latent_vapor_remaining: float):
        """
        Determine the current thermodynamic phase.

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

        Returns
        -------
        str
            Phase: "solid", "fusion", "liquid", "sublimation", or "vapor"
        """
        fusion_T = props["fusion_T"]
        vapor_T = props.get("vapor_T", fusion_T + 500.0)

        # Calculate vapor pressure at current temperature
        Pvap = self.vapor.vapor_pressure(material_key, max(fusion_T, T_surface))

        # Get triple point pressure for this material
        P_triple = self.triple_point_pressure.get(material_key, 100.0)

        # Calculate fusion fraction consumed (with safe division)
        L_fusion_total = mass * props["L_fusion"]
        if L_fusion_total > 1e-9:
            fusion_fraction_consumed = 1.0 - (latent_fusion_remaining / L_fusion_total)
        else:
            fusion_fraction_consumed = 1.0

        # Phase determination logic
        if T_surface < fusion_T:
            # Below melting point
            if ambient_pressure < P_triple and Pvap > ambient_pressure:
                # Low pressure regime: sublimation possible
                phase = "sublimation"
            else:
                phase = "solid"

        elif T_surface < vapor_T:
            # Between fusion and vaporization temperatures
            if ambient_pressure < P_triple:
                # Below triple point: no stable liquid phase
                # Material sublimates directly
                phase = "sublimation"
            elif latent_fusion_remaining > 1e-6 and fusion_fraction_consumed < 0.9:
                # Still consuming fusion latent heat
                phase = "fusion"
            else:
                # Fusion complete, in liquid phase
                phase = "liquid"

        else:
            # Above vaporization temperature
            if Pvap > ambient_pressure:
                # Vapor pressure exceeds ambient: vaporization
                if latent_vapor_remaining > 1e-6:
                    phase = "vapor"
                else:
                    # Latent heat exhausted but conditions favor vaporization
                    phase = "vapor"
            else:
                # High pressure suppresses vaporization
                phase = "liquid"

        self.detector.update(phase)
        return phase

    def is_sublimating(self, phase: str) -> bool:
        """Check if material is undergoing sublimation."""
        return phase in ("sublimation", "vapor")

    def core_sublimated(self) -> bool:
        """
        Check if the core has been in a vaporizing state for sustained period.

        This indicates the meteoroid is being destroyed by sublimation/vaporization.
        """
        return self.detector.sustained()
