import math


class Atmosphere:
    def __init__(self):
        # Sea level reference
        self.rho0 = 1.225      # kg/m^3
        self.scale_height = 8500.0  # m
        self.g0 = 9.80665
        # Named layers (simplified Earth, by altitude in meters)
        self.layers = [
            ("Exosphere", 700_000, float("inf")),
            ("Thermosphere", 85_000, 700_000),
            ("Mesosphere", 50_000, 85_000),
            ("Stratosphere", 11_000, 50_000),
            ("Troposphere", 250, 11_000),
            ("Surface", 0, 250)
        ]

        self.p0 = 101325.0  # Pa, sea-level pressure
        self.T0 = 288.15  # K, sea-level temperature
        self.L = 0.0065  # K/m, lapse rate (troposphere)
        self.R = 287.058  # J/(kg*K)
        self.g0 = 9.80665

    def get_layer_name(self, altitude_m: float) -> str:
        h = max(0.0, altitude_m)
        for name, h_min, h_max in self.layers:
            if h_min <= h < h_max:
                return name
        # Above defined exosphere threshold: deep space
        if h >= self.layers[0][2]:
            return "Deep space"
        return "Ground"

    def get_density(self, altitude_m: float) -> float:
        h = max(0.0, altitude_m)
        # Simple exponential model
        return self.rho0 * math.exp(-h / self.scale_height)

    def get_temperature(self, altitude_m: float) -> float:
        h = max(0.0, altitude_m)
        # ISA-like, simplified piecewise (Kelvin)
        if h < 11_000.0:
            return 288.15 - 0.0065 * h
        elif h < 20_000.0:
            return 216.65
        elif h < 32_000.0:
            return 216.65 + 0.001 * (h - 20_000.0)
        elif h < 47_000.0:
            return 228.65 + 0.0028 * (h - 32_000.0)
        else:
            # Upper-atmosphere simplification
            return 270.0

    def get_gravity(self, altitude_m: float) -> float:
        # Mild decrease with altitude (Earth radius ~ 6.371e6 m)
        R_earth = 6_371_000.0
        r = R_earth + max(0.0, altitude_m)
        return self.g0 * (R_earth / r) ** 2

    def get_pressure(self, altitude_m: float) -> float:
        h = max(0.0, altitude_m)
        # ISA troposphere to 11 km, then simple exponential above
        if h <= 11_000.0:
            # Barometric formula (non-isothermal)
            return self.p0 * (1.0 - self.L * h / self.T0) ** (self.g0 / (self.R * self.L))
        else:
            # Isothermal above 11 km using temperature at 11 km
            T11 = self.T0 - self.L * 11_000.0
            p11 = self.p0 * (1.0 - self.L * 11_000.0 / self.T0) ** (self.g0 / (self.R * self.L))
            return p11 * math.exp(-(h - 11_000.0) * self.g0 / (self.R * T11))
