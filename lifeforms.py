import math


class Lifeform:
    PRESETS = {
        # Name: (critical_temp_K, description)
        "human": (323.0, "Human baseline tolerance"),
        "extremophile": (373.0, "Thermotolerant microbe"),
        "plant": (318.0, "Generalized plant tissue"),
        "ice_creature": (273.0, "Cold-adapted organism"),
    }

    def __init__(self, name: str, critical_temp_K: float, desc: str = ""):
        self.name = name
        self.critical_temp = float(critical_temp_K)
        self.desc = desc

    @classmethod
    def from_preset(cls, key: str):
        k = key.lower()
        if k not in cls.PRESETS:
            raise ValueError(f"Unknown lifeform '{key}'. Choose from: {', '.join(cls.PRESETS.keys())}")
        Tcrit, desc = cls.PRESETS[k]
        return cls(name=k, critical_temp_K=Tcrit, desc=desc)

    def survival_probability(self, temperature_profile, dt: float) -> float:
        """
        Converts thermal exposure above critical temp into a probability in [0,1].
        - Compute exposure integral of (T - Tcrit)^p over time with p>1.
        - Map exposure to survival with an exponential hazard model.
        """
        if not temperature_profile:
            return 1.0

        p = 1.5
        scale = 100.0  # K scaling inside the power
        hazard_rate = 2.5e-3  # tune for reasonable ranges

        exposure = 0.0
        for T in temperature_profile:
            delta = max(0.0, T - self.critical_temp)
            exposure += ((delta / scale) ** p) * dt

        # Survival probability
        prob = math.exp(-hazard_rate * exposure)
        # Clamp to [0,1]
        return max(0.0, min(1.0, prob))