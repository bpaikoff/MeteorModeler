class AblationModel:
    """
    Computes ablation (mass loss) rate based on net heating power and material properties.
    This sits downstream of HeatingModel and PhaseModel.
    """

    def compute_mass_loss(self, P_net: float, props: dict, T_surface: float) -> float:
        """
        Parameters
        ----------
        P_net : float
            Net heating power incident on the body [W].
        props : dict
            Material properties (must include 'L_vapor' or 'L', 'fusion_T', 'vapor_T').
        T_surface : float
            Current surface temperature [K].

        Returns
        -------
        mdot : float
            Mass loss rate [kg/s].
        """

        # If below fusion temperature, no ablation
        if T_surface < 0.9 * props["fusion_T"]:
            return 0.0

        # If between fusion and vaporization thresholds, modest ablation
        if props.get("vapor_T") and T_surface < props["vapor_T"]:
            superheat = max(0.0, T_surface - props["fusion_T"])
            f_ablate = 0.15 + 0.25 * (superheat / (superheat + 500.0))
        else:
            # Above vaporization threshold, ablation dominates
            superheat = max(0.0, T_surface - props.get("vapor_T", props["fusion_T"]))
            f_ablate = 0.5 + 0.35 * (superheat / (superheat + 1000.0))

        f_ablate = min(0.9, f_ablate)

        # Latent heat of vaporization (fall back to L if not split)
        L_vap = props.get("L_vapor", props.get("L", 1e6))

        mdot = (f_ablate * max(0.0, P_net)) / max(1e-9, L_vap)
        return mdot