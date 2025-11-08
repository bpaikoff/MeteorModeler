# energy_balance.py
from dataclasses import dataclass

@dataclass
class EnergyPartition:
    E_sensible: float = 0.0
    E_fusion: float = 0.0
    E_vapor: float = 0.0


class EnergyBalanceModel:
    def partition(self, P_net: float, dt: float, mass: float, T_surface: float, props: dict,
                  ambient_pressure: float, vapor_model, material_key: str,
                  latent_fusion_remaining: float, latent_vapor_remaining: float,
                  min_flux_for_phase: float = 100.0):
        E_avail = max(0.0, P_net) * dt
        cp = props["cp"]
        fusion_T = props["fusion_T"]
        part = EnergyPartition()

        # Sensible heating to fusion
        if T_surface < fusion_T and E_avail > 0.0:
            dT_needed = fusion_T - T_surface
            E_needed = dT_needed * mass * cp
            use = min(E_avail, E_needed)
            part.E_sensible = use
            E_avail -= use

        # Fusion latent (consume remaining inventory)
        if E_avail > 0.0 and latent_fusion_remaining > 0.0:
            tick_cap = 0.35  # fraction cap per tick
            use = min(E_avail, tick_cap * latent_fusion_remaining)
            part.E_fusion = use
            E_avail -= use

        # Vaporization only if sustained thermodynamic drive AND we crossed fusion fully
        Pvap = vapor_model.vapor_pressure(material_key, max(fusion_T, T_surface))
        if (E_avail > 0.0 and latent_fusion_remaining <= 1e-6 and
            P_net >= min_flux_for_phase and Pvap > ambient_pressure and
            latent_vapor_remaining > 0.0):
            use = min(E_avail, latent_vapor_remaining)
            part.E_vapor = use
            E_avail -= use

        return part

    def apply(self, part: EnergyPartition, mass: float, T_surface: float, props: dict):
        cp = props["cp"]
        dT_sensible = part.E_sensible / max(1e-9, mass * cp)
        T_new = T_surface + dT_sensible

        L_vap = props.get("L_vapor", props.get("L", 1e6))
        mass_loss = part.E_vapor / max(1e-9, L_vap)

        return T_new, mass_loss