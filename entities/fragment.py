# fragment.py
class FragmentModel:
    def __init__(self, heating, phase, energy):
        self.heat = heating
        self.phase = phase
        self.energy = energy

    def update_particles(self, meteor, rho, T_air, v, dt, props, ambient_pressure, material_key, angle_rad, latent_fusion_remaining, latent_vapor_remaining):
        if not meteor.particles:
            return

        alive = [p for p in meteor.particles if p.alive]
        for p in alive:
            q_net = self.heat.net_heat_flux(meteor.eq_altitude, rho, v, p.emissivity, p.T, T_air, material_key, angle_rad)
            P_net = q_net * max(1e-6, p.A)

            # Partition energy for this particle
            partitions = self.energy.partition(P_net, dt, p.m, p.T, props, ambient_pressure, self.phase.vapor, material_key, latent_fusion_remaining, latent_vapor_remaining)
            T_new, mass_loss = self.energy.apply(partitions, p.m, p.T, props)

            p.T = T_new
            p.m = max(0.0, p.m - mass_loss)
            p.phase = self.phase.evaluate(p.m, props, p.T, ambient_pressure, material_key, latent_fusion_remaining, latent_vapor_remaining)

            if p.m <= 0.0:
                p.alive = False
                p.r = 0.0
            else:
                p.shrink_from_mass(props["density"])

        # Coupling core to fragment average
        alive = [p for p in meteor.particles if p.alive]
        fragment_view = 0.0
        if alive:
            A_frag = sum(p.A for p in alive)
            fragment_view = min(0.6, A_frag / max(1e-6, meteor.area))
            T_avg_frag = sum(p.T for p in alive) / len(alive)
            T_s, T_c, _ = self.heat.update_surface_and_core_temps(
                altitude_m=meteor.eq_altitude,
                A=max(1e-6, meteor.area), m=max(1e-6, meteor.mass), cp=props["cp"], emissivity=props["emissivity"],
                rho=rho, v=v, T_surface=meteor.temperature, T_core=meteor.core_temperature,
                T_air=T_air, dt=dt,
                composition=material_key,  # or meteor.composition
                angle_rad=angle_rad,  # pass from trajectory
                fragment_view_factor=fragment_view,  # computed from fragment area vs core area
            )

            blend = 0.15
            meteor.temperature = (1 - blend) * T_s + blend * T_avg_frag
            meteor.core_temperature = T_c
            meteor.clamp_temperature()