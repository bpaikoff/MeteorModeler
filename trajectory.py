import math


class Trajectory:
    def __init__(self, initial_altitude: float, initial_velocity: float, entry_angle_deg: float):
        self.altitude = float(initial_altitude)
        self.velocity = max(0.0, float(initial_velocity))
        self.last_terminal_velocity = None
        self.entry_angle_deg = float(entry_angle_deg)
        self.angle_rad = math.radians(self.entry_angle_deg)

    def update(self, meteor, rho, dt: float):
        # Update meteor altitude reference for heating model
        meteor.eq_altitude = self.altitude

        # Areas and drag
        A_core = max(1e-6, meteor.area)
        A_frag = sum(p.A for p in meteor.particles if p.alive) if meteor.particles else 0.0
        A_total = max(1e-6, A_core + A_frag)

        Cd_core = meteor.cd
        Cd_frag = meteor.props["Cd"] * 1.1 if meteor.particles else 0.0

        v = max(0.0, self.velocity)

        # Drag forces
        D_core = 0.5 * rho * v * v * Cd_core * A_core
        D_frag = 0.5 * rho * v * v * Cd_frag * A_frag
        D_total = D_core + D_frag

        # Mass
        m_frag = sum(p.m for p in meteor.particles if p.alive) if meteor.particles else 0.0
        m_total = max(1e-6, meteor.mass + m_frag)

        # Gravity accelerates downward along the path; drag decelerates
        g = 9.81
        a_grav = g * math.sin(self.angle_rad)          # positive along descent
        a_drag = D_total / m_total                     # positive opposing motion

        # Semi-implicit update for stability
        dv_dt = a_grav - a_drag
        v_next = self.velocity + dv_dt * dt

        # Terminal velocity estimate; prevents spurious zeroing aloft
        # v_t = sqrt( (2 m g sin(angle)) / (rho Cd A) ), use total area and effective Cd
        Cd_eff = (Cd_core * A_core + Cd_frag * A_frag) / A_total
        v_term = math.sqrt(max(0.0, (2.0 * m_total * g * math.sin(self.angle_rad)) /
                               max(1e-6, rho * Cd_eff * A_total)))
        self.last_terminal_velocity = v_term

        # Constrain only if drag dominates (a_drag >> a_grav) and v_next would flip sign
        if v_next < 0.0 and a_drag > 2.0 * a_grav:
            v_next = max(0.0, 0.95 * v_term)

        self.velocity = max(0.0, v_next)

        # Altitude update (use updated velocity)
        self.altitude = max(0.0, self.altitude - self.velocity * math.sin(self.angle_rad) * dt)

    def has_hit_ground(self) -> bool:
        return self.altitude <= 0.0