# equations.py
from entities import Compositions, FragmentModel
from physics import (
    HeatingModel,
    VaporPressureModel,
    PhaseModel,
    EnergyBalanceModel,
    AblationModel,
    BreakupModel,
    VaporizationDetector,
    PlasmaModel,
    PlasmaStats,
    HeatFluxBreakdown,
)


class Equations:
    def __init__(self):
        self.materials = Compositions()
        self.heat = HeatingModel()
        self.vapor = VaporPressureModel()
        self.detector = VaporizationDetector(window=50)  # ~1s history at dt=0.02
        self.phase = PhaseModel(self.vapor, self.detector)
        self.energy = EnergyBalanceModel()
        self.ablation = AblationModel()
        self.breakup = BreakupModel()
        self.fragments = FragmentModel(self.heat, self.phase, self.energy)
        self.plasma = PlasmaModel()

        # Latest computed stats for UI access
        self.latest_heat_flux: HeatFluxBreakdown | None = None
        self.latest_plasma_stats: PlasmaStats | None = None

    def get_composition(self, name: str):
        return self.materials.get(name)

    def update_thermo_ablation(self, meteor, rho, T_air, v, dt, ambient_pressure, angle_rad):
        props = meteor.props
        A = max(1e-6, meteor.area)
        m = max(1e-6, meteor.mass)
        material_key = meteor.composition

        # Compute heat flux breakdown for UI
        self.latest_heat_flux = self.heat.compute_heat_flux_breakdown(
            altitude_m=meteor.eq_altitude,
            rho=rho, v=v,
            emissivity_base=props["emissivity"],
            T_surface=meteor.temperature,
            T_air=T_air,
            composition=material_key,
            angle_rad=angle_rad,
            fragment_view_factor=0.0 if not meteor.fragmented else 0.15,
        )

        # Compute plasma stats for UI
        self.latest_plasma_stats = self.plasma.compute_all(
            velocity=v,
            altitude_m=meteor.eq_altitude,
            density=rho,
            temperature_freestream=T_air,
            meteor_radius=meteor.radius,
        )

        # Pre-breakup
        if not meteor.fragmented or not meteor.particles:
            T_s, T_c, P_net = self.heat.update_surface_and_core_temps(
                altitude_m=meteor.eq_altitude,
                A=A, m=m, cp=props["cp"], emissivity=props["emissivity"],
                rho=rho, v=v, T_surface=meteor.temperature, T_core=meteor.core_temperature,
                T_air=T_air, dt=dt,
                composition=meteor.composition,
                angle_rad=angle_rad,
                radius=meteor.radius,
                density=props["density"],
            )

            # Partition energy
            part = self.energy.partition(
                P_net, dt, m, T_s, props, ambient_pressure, self.vapor, material_key,
                meteor.latent_fusion_remaining, meteor.latent_vapor_remaining
            )
            T_s_new, mass_loss_energy = self.energy.apply(part, m, T_s, props)

            # Update latent inventories
            meteor.latent_fusion_remaining = max(0.0, meteor.latent_fusion_remaining - part.E_fusion)
            meteor.cum_fusion_used += part.E_fusion
            meteor.latent_vapor_remaining = max(0.0, meteor.latent_vapor_remaining - part.E_vapor)
            meteor.cum_vapor_used += part.E_vapor

            # Phase evaluation (pass P_net for energy-limited sublimation)
            phase_now = self.phase.evaluate(
                m, props, T_s_new, ambient_pressure, material_key,
                meteor.latent_fusion_remaining, meteor.latent_vapor_remaining,
                P_net=P_net
            )
            meteor.phase = phase_now
            meteor.core_phase = phase_now
            meteor.core_sublimated = self.phase.core_sublimated()

            # Ablation
            mdot = self.ablation.compute_mass_loss(P_net, props, T_s_new)

            # Update meteor
            meteor.temperature = T_s_new
            meteor.core_temperature = T_c
            meteor.mass = max(0.0, meteor.mass - (mass_loss_energy + mdot * dt))
            meteor.update_geometry_from_mass()
            meteor.clamp_temperature()
            return

        # Post-breakup
        self.fragments.update_particles(meteor, rho, T_air, v, dt, props,
                                        ambient_pressure, material_key, angle_rad, meteor.latent_fusion_remaining, meteor.latent_vapor_remaining)

    def check_breakup(self, meteor, rho, v, angle_rad):
        if self.breakup.should_breakup(meteor, rho, v, angle_rad):
            # Use progressive fragmentation with computed fragment count
            n_fragments = self.breakup.compute_fragment_count(meteor, rho, v, angle_rad)
            if n_fragments > 0:
                self.breakup.fragment_cascade(meteor, n_fragments)