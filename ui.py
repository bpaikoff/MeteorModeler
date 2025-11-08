import pygame
import sys
import math


class UI:
    def __init__(self, width=900, height=700):
        self.width = width
        self.height = height
        self.screen = None
        self.clock = None
        self.font = None
        self.large_font = None
        self.max_alt_km = 100.0

        # Camera and scale
        self.margin_top = 40
        self.margin_bottom = 60
        self.sky_color = (10, 16, 32)
        self.ground_color = (45, 80, 45)

        # Trail
        self.trail = []

    def display_intro(self):
        print("=== Meteor Descent Survival Simulator ===")
        print("You’ll pick a meteor composition and diameter, and a lifeform.")
        print("The simulation will show a falling meteor; after it ends,")
        print("you’ll see a survival chance based on thermal exposure.\n")

    def get_user_choices(self):
        # Composition
        valid = ["iron", "stone", "ice", "carbonaceous"]
        comp = input(f"Choose composition {valid} [stone]: ").strip().lower() or "stone"
        while comp not in valid:
            comp = input(f"Please choose one of {valid}: ").strip().lower()

        # Diameter
        while True:
            s = input("Enter meteor diameter in meters (e.g., 1 for 1 m, 10, 50) [5]: ").strip()
            if s == "":
                diameter = 5.0
                break
            try:
                diameter = float(s)
                if diameter <= 0.0:
                    raise ValueError
                break
            except Exception:
                print("Please enter a positive number.")

        # Lifeform
        forms = ["human", "extremophile", "plant", "ice_creature"]
        life = input(f"Choose lifeform {forms} [extremophile]: ").strip().lower() or "extremophile"
        while life not in forms:
            life = input(f"Please choose one of {forms}: ").strip().lower()

        print(f"\nSelection → Composition: {comp}, Diameter: {diameter} m, Lifeform: {life}\n")
        return comp, diameter, life

    def init_scene(self, max_alt_km=100.0):
        pygame.init()
        pygame.display.set_caption("Meteor Descent Survival Simulator")
        self.screen = pygame.display.set_mode((self.width, self.height))
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont("consolas", 18)
        self.large_font = pygame.font.SysFont("consolas", 28, bold=True)
        self.max_alt_km = max_alt_km
        self.trail.clear()

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False
            if event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE:
                return False
        return True

    def altitude_to_screen_y(self, altitude_m):
        alt_km = max(0.0, altitude_m) / 1000.0
        # Map [0 km, max_alt_km] to [ground_y, margin_top]
        ground_y = self.height - self.margin_bottom
        top_y = self.margin_top
        t = min(1.0, alt_km / self.max_alt_km)
        y = ground_y - (ground_y - top_y) * t
        return int(y)

    def temp_to_color(self, T):
        # Map 300K..3000K to blue→red→white
        T0, T1 = 300.0, 3000.0
        x = (max(T0, min(T, T1)) - T0) / (T1 - T0)
        # Simple gradient
        r = int(100 + 155 * x)
        g = int(80 + 120 * (1 - abs(2 * x - 1)))
        b = int(120 + 80 * (1 - x))
        return (r, g, b)

    def update_display(self, meteor, traj, atm, time_elapsed=0.0):
        self.screen.fill(self.sky_color)

        pygame.draw.rect(
            self.screen,
            self.ground_color,
            pygame.Rect(0, self.height - self.margin_bottom, self.width, self.margin_bottom),
        )

        y = self.altitude_to_screen_y(traj.altitude)
        x = int(self.width * 0.5)

        self.trail.append((x, y))
        if len(self.trail) > 200:
            self.trail.pop(0)
        for i, (tx, ty) in enumerate(self.trail):
            alpha = i / len(self.trail)
            color = (200, int(150 * alpha), int(80 * alpha))
            pygame.draw.circle(self.screen, color, (tx, ty), 2)

        # Outer-temperature color for visuals
        color = self.temp_to_color(meteor.temperature)
        radius_px = max(2, int(5 + meteor.radius * 0.5))
        pygame.draw.circle(self.screen, color, (x, y), radius_px)

        # Animate fragments
        if meteor.particles:
            for p in meteor.particles:
                if not p.alive:
                    continue
                px = x + int(p.x * 0.02)
                py = y + int(p.x * 0.005)
                p_color = self.temp_to_color(p.T)
                p_radius = max(1, int(2 + p.r * 0.3))
                pygame.draw.circle(self.screen, p_color, (px, py), p_radius)

        # HUD: conversions, particle count, time, layer
        def K_to_F(Tk):
            return (Tk - 273.15) * 9 / 5 + 32

        particle_count = 0
        if meteor.particles:
            particle_count = max(1, sum(1 for p in meteor.particles if p.alive))

        layer_name = atm.get_layer_name(traj.altitude)
        v_term = traj.last_terminal_velocity or 0.0
        core_status = "sublimated" if getattr(meteor, "core_sublimated", False) else getattr(meteor, "core_phase",
                                                                                             "n/a")
        hud_lines = [
            f"Time: {time_elapsed:6.2f} s",
            f"Layer: {layer_name}",
            f"Alt: {traj.altitude / 1000:.1f} km",
            f"Vel: {traj.velocity:.0f} m/s",
            f"Term Vel: {v_term:.0f} m/s",
            f"Phase: {meteor.phase} | Core: {core_status}",
            f"Outer Temp: {meteor.temperature:.0f} K / {K_to_F(meteor.temperature):.0f} °F",
            f"Core Temp: {meteor.core_temperature:.0f} K / {K_to_F(meteor.core_temperature):.0f} °F",
            f"Mass: {meteor.mass:.1f} kg",
            f"Diam: {meteor.diameter:.2f} m",
            f"Particles: {particle_count}",
        ]

        for i, line in enumerate(hud_lines):
            surf = self.font.render(line, True, (255, 255, 255))
            self.screen.blit(surf, (10, 10 + 20 * i))

    def show_result_overlay(self, text):
        # Semi-transparent overlay
        overlay = pygame.Surface((self.width, self.height), pygame.SRCALPHA)
        overlay.fill((0, 0, 0, 180))
        self.screen.blit(overlay, (0, 0))

        # Centered text
        lines = text.split("\n")
        for i, line in enumerate(lines):
            surf = self.large_font.render(line, True, (255, 255, 255))
            rect = surf.get_rect(center=(self.width // 2, self.height // 2 + i * 40))
            self.screen.blit(surf, rect)

    def flip(self):
        pygame.display.flip()
        self.clock.tick(60)

    def shutdown(self):
        pygame.quit()
        sys.exit()