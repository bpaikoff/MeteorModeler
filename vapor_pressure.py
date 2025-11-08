# vapor_pressure.py
import math

class VaporPressureModel:
    def __init__(self):
        # Latent heats (J/kg)
        self.L_sublimation = 2.83e6   # ice → vapor
        self.L_vaporization = 2.26e6  # liquid → vapor
        self.Rv = 461.5               # J/(kg*K), water vapor gas constant

        # Reference triple point
        self.T0 = 273.16  # K
        self.P0 = 611.0   # Pa

    def vapor_pressure(self, material: str, T_K: float) -> float:
        """
        Returns equilibrium vapor pressure [Pa] at temperature T_K.
        Uses Clausius–Clapeyron anchored at triple point.
        """
        if material != "ice":
            # fallback for non-volatiles
            if T_K < 1000.0:
                return 1.0
            return 10.0 ** (0.01 * (T_K - 1000.0))

        # Choose latent heat depending on phase range
        if T_K < self.T0:
            L = self.L_sublimation
        else:
            L = self.L_vaporization

        exponent = -L / self.Rv * (1.0 / T_K - 1.0 / self.T0)
        P = self.P0 * math.exp(exponent)
        return max(1e-6, P)  # clamp to avoid zero