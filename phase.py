# phase.py
class PhaseModel:
    def __init__(self, vapor_model, detector):
        self.vapor = vapor_model
        self.detector = detector

    def evaluate(self, mass, props, T_surface, ambient_pressure, material_key: str,
                 latent_fusion_remaining: float, latent_vapor_remaining: float):
        fusion_T = props["fusion_T"]
        Pvap = self.vapor.vapor_pressure(material_key, max(fusion_T, T_surface))

        if T_surface < fusion_T:
            phase = "solid"
        elif  latent_fusion_remaining / (mass * props["L_fusion"]) < 0.1:
            # Allow transition to liquid once a significant fraction of fusion latent is consumed
            phase = "liquid"
        elif latent_fusion_remaining > 1e-6:
            phase = "fusion"  # actively consuming fusion latent
        elif ambient_pressure >= Pvap:
            phase = "liquid"
        else:
            # Only consider vapor phase if there's vapor latent to consume
            phase = "vapor" if latent_vapor_remaining > 1e-6 else "liquid"

        self.detector.update(phase)
        return phase

    def core_sublimated(self) -> bool:
        return self.detector.sustained()