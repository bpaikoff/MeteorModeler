import math
from entities.particles import Particle


class Meteor:
    def __init__(self, composition: str, diameter: float, equations):
        """
        Represents a meteoroid body entering the atmosphere.
        - composition: string key (iron, stone, ice, carbonaceous)
        - diameter: initial diameter in meters
        - equations: Equations coordinator (for composition lookup)
        """
        self.composition = composition.lower()
        self.eq = equations
        self.props = self.eq.get_composition(self.composition)

        self.diameter = max(0.01, float(diameter))
        self.radius = 0.5 * self.diameter
        self.volume = (4.0 / 3.0) * math.pi * self.radius ** 3
        self.mass = self.props["density"] * self.volume

        # Temperatures
        self.temperature = 300.0       # outer/surface K
        self.core_temperature = 300.0  # inner/core K

        # Phase tracking
        self.phase = "solid"
        self.core_phase = "solid"
        self.core_sublimated = False   # triggers sim end if True
        self.latent_fusion_remaining = self.mass * self.props["L_fusion"]  # J
        self.latent_vapor_remaining = self.mass * self.props["L_vapor"]  # J
        self.cum_fusion_used = 0.0
        self.cum_vapor_used = 0.0

        # Fragmentation
        self.fragmented = False
        self.particles = None          # list[Particle] after breakup
        self.area_multiplier = 1.0
        self.cd_multiplier = 1.0

        # Rendering helpers
        self.max_temp_seen = self.temperature
        self.burned_out = False

        # Altitude reference for heating model (set by trajectory each tick)
        self.eq_altitude = 120_000.0

    @property
    def area(self):
        # Effective cross-sectional area
        return math.pi * (self.radius ** 2) * self.area_multiplier

    @property
    def cd(self):
        return self.props["Cd"] * self.cd_multiplier

    def update_geometry_from_mass(self):
        self.mass = max(self.mass, 0.0)
        if self.mass <= 0.0:
            self.radius = 0.0
            self.diameter = 0.0
            self.volume = 0.0
            self.burned_out = True
            return

        self.volume = self.mass / self.props["density"]
        self.radius = ((3.0 * self.volume) / (4.0 * math.pi)) ** (1.0 / 3.0)
        self.diameter = 2.0 * self.radius

    def note_fragmentation(self):
        if not self.fragmented:
            self.fragmented = True
            # Convert bulk body into N particles
            N = 64
            density = self.props["density"]
            total_mass = self.mass
            base_r = max(0.02, 0.15 * self.radius)
            masses = [total_mass / N for _ in range(N)]
            self.particles = []
            for i in range(N):
                xi = (i / (N - 1))
                x_pos = -3.0 * self.radius * xi
                p = Particle(
                    x=x_pos,
                    r=base_r,
                    m=masses[i],
                    T=self.temperature,
                    Cd=self.props["Cd"] * 1.1,
                    emissivity=self.props["emissivity"],
                    cp=self.props["cp"],
                    L=self.props["L_vapor"] if "L_vapor" in self.props else self.props["L"],
                )
                p.shrink_from_mass(density)
                p.phase = "solid"
                self.particles.append(p)

            # Reduce bulk body to a core fraction
            core_fraction = 0.12
            self.mass = total_mass * core_fraction
            self.update_geometry_from_mass()
            self.area_multiplier *= 1.6
            self.cd_multiplier *= 1.15

    def effective_area(self):
        if self.particles:
            return self.area + sum(p.A for p in self.particles if p.alive)
        return self.area

    def is_fully_ablated(self):
        if self.burned_out or self.mass <= 0.0 or self.radius <= 0.0:
            return True
        if self.particles:
            alive_mass = sum(p.m for p in self.particles if p.alive)
            return alive_mass <= 0.0 and self.mass <= 0.0
        return False

    def clamp_temperature(self):
        self.temperature = max(100.0, min(self.temperature, 6000.0))
        self.core_temperature = max(100.0, min(self.core_temperature, 6000.0))
        self.max_temp_seen = max(self.max_temp_seen, self.temperature)