import pygame
import sys
import math


class UI:
    def __init__(self, width=1280, height=700):
        self.width = width
        self.height = height
        self.screen = None
        self.clock = None
        self.font = None
        self.small_font = None
        self.large_font = None
        self.max_alt_km = 100.0

        # Camera and scale
        self.margin_top = 40
        self.margin_bottom = 60
        self.sky_color = (10, 16, 32)
        self.ground_color = (45, 80, 45)

        # Layout: left panel (200px) | visualization | right panel (320px)
        self.left_panel_width = 200
        self.right_panel_width = 320
        self.viz_left = self.left_panel_width
        self.viz_right = self.width - self.right_panel_width

        # Trail
        self.trail = []

        # Pause state
        self.paused = False

        # Colors for panels
        self.panel_bg = (20, 25, 40)
        self.panel_border = (60, 70, 100)
        self.text_color = (220, 220, 220)
        self.header_color = (100, 180, 255)
        self.value_color = (180, 255, 180)
        self.bar_bg = (40, 45, 60)
        self.bar_conv = (255, 150, 80)    # Orange for convective
        self.bar_rad = (150, 100, 255)    # Purple for radiative
        self.bar_cat = (80, 200, 150)     # Teal for catalytic
        self.bar_cool = (100, 150, 255)   # Blue for cooling

    def display_intro(self):
        print("=== Meteor Descent Survival Simulator ===")
        print("You'll pick a meteor composition and diameter, and a lifeform.")
        print("The simulation will show a falling meteor; after it ends,")
        print("you'll see a survival chance based on thermal exposure.\n")

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
        self.font = pygame.font.SysFont("consolas", 16)
        self.small_font = pygame.font.SysFont("consolas", 13)
        self.large_font = pygame.font.SysFont("consolas", 28, bold=True)
        self.max_alt_km = max_alt_km
        self.trail.clear()

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    return False
                if event.key == pygame.K_SPACE or event.key == pygame.K_p:
                    self.paused = not self.paused
        return True

    def is_paused(self):
        return self.paused

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

    def draw_panel_background(self, x, y, width, height, title=None):
        """Draw a panel with background and optional title."""
        # Background
        pygame.draw.rect(self.screen, self.panel_bg, pygame.Rect(x, y, width, height))
        # Border
        pygame.draw.rect(self.screen, self.panel_border, pygame.Rect(x, y, width, height), 1)
        # Title
        if title:
            surf = self.font.render(title, True, self.header_color)
            self.screen.blit(surf, (x + 8, y + 5))
            # Title underline
            pygame.draw.line(self.screen, self.panel_border,
                           (x + 5, y + 24), (x + width - 5, y + 24), 1)
            return y + 30  # Return y position after title
        return y + 5

    def format_scientific(self, value, precision=2):
        """Format a number in scientific notation for display."""
        if value == 0 or value == float('inf') or math.isnan(value):
            return "0"
        if abs(value) < 1e-20:
            return "~0"
        exp = int(math.floor(math.log10(abs(value))))
        mantissa = value / (10 ** exp)
        if exp == 0:
            return f"{value:.{precision}f}"
        elif abs(exp) <= 3:
            return f"{value:.{precision}f}"
        else:
            return f"{mantissa:.{precision}f}e{exp}"

    def draw_percentage_bar(self, x, y, width, percentages, colors, height=12):
        """Draw a stacked percentage bar."""
        # Background
        pygame.draw.rect(self.screen, self.bar_bg, pygame.Rect(x, y, width, height))

        # Draw segments
        current_x = x
        for pct, color in zip(percentages, colors):
            seg_width = int(width * pct / 100.0)
            if seg_width > 0:
                pygame.draw.rect(self.screen, color,
                               pygame.Rect(current_x, y, seg_width, height))
                current_x += seg_width

        # Border
        pygame.draw.rect(self.screen, self.panel_border, pygame.Rect(x, y, width, height), 1)

    def draw_heat_flux_panel(self, heat_flux, x, y, width):
        """Draw the heat flux breakdown panel."""
        panel_height = 180
        content_y = self.draw_panel_background(x, y, width, panel_height, "Heat Flux Breakdown")

        if heat_flux is None:
            surf = self.small_font.render("No data", True, self.text_color)
            self.screen.blit(surf, (x + 10, content_y + 5))
            return y + panel_height + 5

        line_height = 16

        # Flow regime
        regime_color = (100, 200, 100) if heat_flux.flow_regime == "continuum" else (200, 150, 100)
        surf = self.small_font.render(f"Flow: {heat_flux.flow_regime}", True, regime_color)
        self.screen.blit(surf, (x + 10, content_y))
        content_y += line_height

        # Heat flux values
        def draw_flux_line(label, value, unit="MW/m²"):
            nonlocal content_y
            val_str = self.format_scientific(value / 1e6, 2)  # Convert W/m² to MW/m²
            surf = self.small_font.render(f"{label}:", True, self.text_color)
            self.screen.blit(surf, (x + 10, content_y))
            surf = self.small_font.render(f"{val_str} {unit}", True, self.value_color)
            self.screen.blit(surf, (x + 130, content_y))
            content_y += line_height

        draw_flux_line("Convective", heat_flux.q_convective)
        draw_flux_line("Radiative", heat_flux.q_radiative_shock)
        draw_flux_line("Catalytic", heat_flux.q_catalytic)
        draw_flux_line("Cooling", heat_flux.q_radiative_cooling)
        draw_flux_line("Net", heat_flux.q_net)

        content_y += 5

        # Percentage breakdown label
        surf = self.small_font.render("Contribution breakdown:", True, self.text_color)
        self.screen.blit(surf, (x + 10, content_y))
        content_y += line_height

        # Stacked bar
        self.draw_percentage_bar(
            x + 10, content_y, width - 20,
            [heat_flux.pct_convective, heat_flux.pct_radiative, heat_flux.pct_catalytic],
            [self.bar_conv, self.bar_rad, self.bar_cat]
        )
        content_y += 18

        # Legend
        legend_items = [
            ("Conv", self.bar_conv, heat_flux.pct_convective),
            ("Rad", self.bar_rad, heat_flux.pct_radiative),
            ("Cat", self.bar_cat, heat_flux.pct_catalytic),
        ]
        legend_x = x + 10
        for label, color, pct in legend_items:
            pygame.draw.rect(self.screen, color, pygame.Rect(legend_x, content_y, 10, 10))
            surf = self.small_font.render(f"{label}:{pct:.0f}%", True, self.text_color)
            self.screen.blit(surf, (legend_x + 14, content_y - 2))
            legend_x += 80

        return y + panel_height + 5

    def draw_plasma_panel(self, plasma, x, y, width):
        """Draw the plasma properties panel."""
        panel_height = 310
        content_y = self.draw_panel_background(x, y, width, panel_height, "Plasma Properties")

        if plasma is None:
            surf = self.small_font.render("No data", True, self.text_color)
            self.screen.blit(surf, (x + 10, content_y + 5))
            return y + panel_height + 5

        line_height = 15

        def draw_stat(label, value, unit=""):
            nonlocal content_y
            if isinstance(value, float):
                if value == float('inf') or math.isnan(value):
                    val_str = "N/A"
                elif abs(value) > 1e6 or (abs(value) < 0.01 and value != 0):
                    val_str = self.format_scientific(value, 2)
                else:
                    val_str = f"{value:.2f}"
            else:
                val_str = str(value)
            surf = self.small_font.render(f"{label}:", True, self.text_color)
            self.screen.blit(surf, (x + 10, content_y))
            surf = self.small_font.render(f"{val_str} {unit}", True, self.value_color)
            self.screen.blit(surf, (x + 150, content_y))
            content_y += line_height

        def draw_section(title):
            nonlocal content_y
            content_y += 3
            surf = self.small_font.render(title, True, self.header_color)
            self.screen.blit(surf, (x + 10, content_y))
            content_y += line_height

        # Shock section
        draw_section("-- Shock Layer --")
        draw_stat("Mach number", plasma.mach_number)
        draw_stat("Shock temp", plasma.shock_temperature, "K")
        draw_stat("Density ratio", plasma.shock_density_ratio)
        draw_stat("Standoff dist", plasma.shock_standoff * 100, "cm")

        # Ionization section
        draw_section("-- Ionization --")
        draw_stat("Ion. fraction", plasma.ionization_fraction * 100, "%")
        draw_stat("Electron n_e", plasma.electron_density, "m⁻³")
        draw_stat("Dissoc. O₂", plasma.dissociation_O2 * 100, "%")
        draw_stat("Dissoc. N₂", plasma.dissociation_N2 * 100, "%")

        # Plasma parameters section
        draw_section("-- Plasma Parameters --")
        draw_stat("Debye length", plasma.debye_length * 1e6, "μm")
        draw_stat("Plasma freq", plasma.plasma_frequency / 1e9, "GHz")
        draw_stat("Elec. cond.", plasma.electrical_conductivity, "S/m")
        draw_stat("Mag. Reynolds", plasma.magnetic_reynolds)

        # Flow section
        draw_section("-- Flow --")
        draw_stat("Knudsen", plasma.knudsen_number)
        draw_stat("Reynolds", plasma.reynolds_number)

        return y + panel_height + 5

    def update_display(self, meteor, traj, atm, time_elapsed=0.0, heat_flux=None, plasma_stats=None, breakup_events=0):
        self.screen.fill(self.sky_color)

        # Draw visualization area background
        viz_width = self.viz_right - self.viz_left
        pygame.draw.rect(
            self.screen,
            self.ground_color,
            pygame.Rect(self.viz_left, self.height - self.margin_bottom,
                       viz_width, self.margin_bottom),
        )

        # Meteor position in visualization area
        y = self.altitude_to_screen_y(traj.altitude)
        x = int(self.viz_left + viz_width * 0.5)

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

        # Left panel: HUD
        def K_to_F(Tk):
            return (Tk - 273.15) * 9 / 5 + 32

        # Particle count: 1 (the meteor itself) + any fragments
        if meteor.particles:
            fragment_count = sum(1 for p in meteor.particles if p.alive)
            particle_count = 1 + fragment_count  # Core + fragments
        else:
            particle_count = 1  # Just the meteor body

        layer_name = atm.get_layer_name(traj.altitude)
        v_term = traj.last_terminal_velocity or 0.0
        core_status = "sublimated" if getattr(meteor, "core_sublimated", False) else getattr(meteor, "core_phase",
                                                                                             "n/a")

        # Draw left panel background
        panel_y = self.draw_panel_background(0, 0, self.left_panel_width, 300, "Meteor Status")

        hud_lines = [
            f"Time: {time_elapsed:6.2f} s",
            f"Layer: {layer_name}",
            f"Alt: {traj.altitude / 1000:.1f} km",
            f"Vel: {traj.velocity:.0f} m/s",
            f"Term Vel: {v_term:.0f} m/s",
            f"Phase: {meteor.phase}",
            f"Core: {core_status}",
            f"Surf: {meteor.temperature:.0f} K",
            f"      ({K_to_F(meteor.temperature):.0f} °F)",
            f"Core: {meteor.core_temperature:.0f} K",
            f"      ({K_to_F(meteor.core_temperature):.0f} °F)",
            f"Mass: {meteor.mass:.1f} kg",
            f"Diam: {meteor.diameter:.2f} m",
            f"Particles: {particle_count}",
            f"Breakups: {breakup_events}",
        ]

        for i, line in enumerate(hud_lines):
            surf = self.small_font.render(line, True, self.text_color)
            self.screen.blit(surf, (8, panel_y + 15 * i))

        # Right panel: Heat flux and plasma stats
        right_x = self.viz_right + 5
        panel_y = 0

        panel_y = self.draw_heat_flux_panel(heat_flux, right_x, panel_y, self.right_panel_width - 10)
        panel_y = self.draw_plasma_panel(plasma_stats, right_x, panel_y, self.right_panel_width - 10)

        # Draw pause indicator if paused
        if self.paused:
            self.draw_pause_indicator()

        # Draw controls hint at bottom of visualization area
        self.draw_controls_hint()

    def draw_pause_indicator(self):
        """Draw pause indicator in center of visualization area."""
        viz_center_x = (self.viz_left + self.viz_right) // 2
        viz_center_y = self.height // 2

        # Semi-transparent background
        overlay_width = 200
        overlay_height = 80
        overlay = pygame.Surface((overlay_width, overlay_height), pygame.SRCALPHA)
        overlay.fill((0, 0, 0, 160))
        self.screen.blit(overlay, (viz_center_x - overlay_width // 2,
                                   viz_center_y - overlay_height // 2))

        # Pause icon (two vertical bars)
        bar_width = 12
        bar_height = 40
        bar_gap = 16
        bar_color = (255, 255, 255)
        left_bar_x = viz_center_x - bar_gap // 2 - bar_width
        right_bar_x = viz_center_x + bar_gap // 2
        bar_y = viz_center_y - bar_height // 2 - 10

        pygame.draw.rect(self.screen, bar_color,
                        pygame.Rect(left_bar_x, bar_y, bar_width, bar_height))
        pygame.draw.rect(self.screen, bar_color,
                        pygame.Rect(right_bar_x, bar_y, bar_width, bar_height))

        # PAUSED text
        surf = self.font.render("PAUSED", True, (255, 255, 255))
        rect = surf.get_rect(center=(viz_center_x, viz_center_y + 25))
        self.screen.blit(surf, rect)

    def draw_controls_hint(self):
        """Draw controls hint at bottom of screen."""
        hint_text = "[Space/P] Pause  |  [Esc] Quit"
        surf = self.small_font.render(hint_text, True, (150, 150, 150))
        self.screen.blit(surf, (self.viz_left + 10, self.height - 25))

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
