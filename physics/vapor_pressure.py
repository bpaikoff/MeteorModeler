# vapor_pressure.py
"""
Vapor pressure models for meteoroid materials using Clausius-Clapeyron equation.

The Clausius-Clapeyron equation relates vapor pressure to temperature:
    ln(P/P0) = -(L/R) * (1/T - 1/T0)

where:
    L = latent heat of vaporization (J/mol)
    R = universal gas constant (8.314 J/(mol·K))
    T0, P0 = reference temperature and pressure

For meteor entry, vapor pressure determines when material sublimates or
vaporizes. At low ambient pressures (upper atmosphere), materials sublimate
when their vapor pressure exceeds ambient pressure.
"""
import math


class VaporPressureModel:
    def __init__(self):
        # Universal gas constant [J/(mol·K)]
        self.R = 8.314

        # Material-specific parameters
        # Format: {material: (T_ref, P_ref, L_vap_molar, molar_mass, description)}
        # T_ref [K], P_ref [Pa], L_vap [J/mol], M [kg/mol]
        self.materials = {
            # Water ice: well-characterized triple point
            # L_sublimation = 51.06 kJ/mol, L_vaporization = 40.65 kJ/mol
            "ice": {
                "T_triple": 273.16,  # K
                "P_triple": 611.657,  # Pa
                "L_sublimation": 51060.0,  # J/mol (ice → vapor)
                "L_vaporization": 40650.0,  # J/mol (liquid → vapor)
                "molar_mass": 0.018015,  # kg/mol
            },

            # Iron: vaporization at high temperatures
            # Boiling point ~3134 K at 1 atm, L_vap ≈ 340 kJ/mol
            "iron": {
                "T_ref": 3134.0,  # K (boiling point at 1 atm)
                "P_ref": 101325.0,  # Pa
                "L_vaporization": 340000.0,  # J/mol
                "molar_mass": 0.05585,  # kg/mol
            },

            # Silicate rock (approximated as enstatite MgSiO3)
            # Sublimation/vaporization ~2500-3000 K under vacuum
            "stone": {
                "T_ref": 2800.0,  # K (approximate vaporization)
                "P_ref": 101325.0,  # Pa
                "L_vaporization": 400000.0,  # J/mol (silicate vaporization)
                "molar_mass": 0.1004,  # kg/mol (MgSiO3)
            },

            # Carbonaceous chondrite: complex mixture with volatiles
            # Contains organics that volatilize at lower temps, plus silicates
            # Use weighted-average behavior
            "carbonaceous": {
                "T_ref": 2200.0,  # K (lower due to volatile organics)
                "P_ref": 101325.0,  # Pa
                "L_vaporization": 250000.0,  # J/mol (effective)
                "molar_mass": 0.050,  # kg/mol (effective)
            },
        }

    def vapor_pressure(self, material: str, T_K: float) -> float:
        """
        Returns equilibrium vapor pressure [Pa] at temperature T_K.

        Uses Clausius-Clapeyron equation:
            P = P_ref * exp(-(L/R) * (1/T - 1/T_ref))

        Parameters
        ----------
        material : str
            Material type: 'ice', 'iron', 'stone', 'carbonaceous'
        T_K : float
            Temperature in Kelvin

        Returns
        -------
        float
            Vapor pressure in Pascals
        """
        T_K = max(1.0, T_K)  # Avoid division by zero

        if material not in self.materials:
            # Unknown material: use generic high-temperature vaporization
            return self._generic_vapor_pressure(T_K)

        mat = self.materials[material]

        if material == "ice":
            return self._ice_vapor_pressure(T_K, mat)
        else:
            return self._metal_silicate_vapor_pressure(T_K, mat)

    def _ice_vapor_pressure(self, T_K: float, params: dict) -> float:
        """
        Water ice vapor pressure with sublimation/vaporization transition.

        Below triple point (273.16 K): sublimation (ice → vapor)
        Above triple point: vaporization (liquid → vapor)
        """
        T_triple = params["T_triple"]
        P_triple = params["P_triple"]

        if T_K < T_triple:
            # Sublimation regime
            L = params["L_sublimation"]
        else:
            # Liquid-vapor equilibrium
            L = params["L_vaporization"]

        exponent = -(L / self.R) * (1.0 / T_K - 1.0 / T_triple)
        P = P_triple * math.exp(exponent)

        return max(1e-10, P)

    def _metal_silicate_vapor_pressure(self, T_K: float, params: dict) -> float:
        """
        Vapor pressure for metals and silicates.

        These materials have very low vapor pressure at moderate temperatures
        and only become significant near their boiling points.
        """
        T_ref = params["T_ref"]
        P_ref = params["P_ref"]
        L = params["L_vaporization"]

        # Below ~60% of reference temperature, vapor pressure is negligible
        if T_K < 0.6 * T_ref:
            return 1e-10

        exponent = -(L / self.R) * (1.0 / T_K - 1.0 / T_ref)
        P = P_ref * math.exp(exponent)

        return max(1e-10, min(P, 1e9))  # Clamp to reasonable range

    def _generic_vapor_pressure(self, T_K: float) -> float:
        """
        Generic fallback for unknown materials.
        Assumes high-temperature vaporization similar to silicates.
        """
        if T_K < 1500.0:
            return 1e-10

        # Approximate: P ~ 1 Pa at 2000 K, increasing exponentially
        T_ref = 2500.0
        L_approx = 350000.0  # J/mol
        P_ref = 101325.0

        exponent = -(L_approx / self.R) * (1.0 / T_K - 1.0 / T_ref)
        P = P_ref * math.exp(exponent)

        return max(1e-10, min(P, 1e9))

    def sublimation_temperature(self, material: str, P_ambient: float) -> float:
        """
        Returns the temperature at which vapor pressure equals ambient pressure.
        This is the sublimation/boiling point at the given pressure.

        Parameters
        ----------
        material : str
            Material type
        P_ambient : float
            Ambient pressure [Pa]

        Returns
        -------
        float
            Sublimation/boiling temperature [K]
        """
        P_ambient = max(1e-10, P_ambient)

        if material not in self.materials:
            # Generic: assume ~2500 K at 1 atm, scales with pressure
            return 2500.0 * (1.0 - 0.1 * math.log10(101325.0 / P_ambient))

        mat = self.materials[material]

        if material == "ice":
            T_ref = mat["T_triple"]
            P_ref = mat["P_triple"]
            # Use sublimation latent heat for low pressures
            L = mat["L_sublimation"] if P_ambient < P_ref else mat["L_vaporization"]
        else:
            T_ref = mat["T_ref"]
            P_ref = mat["P_ref"]
            L = mat["L_vaporization"]

        # Invert Clausius-Clapeyron: T = 1 / (1/T_ref - R/L * ln(P/P_ref))
        ln_ratio = math.log(P_ambient / P_ref)
        inv_T = (1.0 / T_ref) - (self.R / L) * ln_ratio

        if inv_T <= 0:
            return 10000.0  # Very high temperature limit

        return 1.0 / inv_T
