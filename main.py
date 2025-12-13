from entities import Meteor
from physics import Atmosphere, Trajectory
from equations import Equations
from lifeforms import Lifeform
from ui import UI
import math


def main():
    ui = UI(width=900, height=700)
    ui.display_intro()

    composition, diameter_m, life_choice = ui.get_user_choices()

    eq = Equations()
    atm = Atmosphere()
    meteor = Meteor(composition=composition, diameter=diameter_m, equations=eq)

    # Entry interface: 120 km, shallow-to-moderate angle, typical NEAs ~ 11–72 km/s; pick 19 km/s for stability
    traj = Trajectory(
        initial_altitude=120_000.0,    # 120 km
        initial_velocity=19_000.0,     # m/s
        entry_angle_deg=45.0
    )

    life = Lifeform.from_preset(life_choice)

    dt = 0.02
    max_sim_time = 900.0  # 15 minutes cap (should finish far sooner)
    t = 0.0

    temperature_profile = []
    ui.init_scene(max_alt_km=120.0)

    running = True
    sim_active = True
    result_text = None
    survival_percent = 0.0

    while running:
        running = ui.handle_events()

        if sim_active:
            rho = atm.get_density(traj.altitude)
            T_air = atm.get_temperature(traj.altitude)
            pressure = atm.get_pressure(traj.altitude)

            eq.update_thermo_ablation(meteor, rho, T_air, traj.velocity, dt, pressure, traj.angle_rad)
            eq.check_breakup(meteor, rho, traj.velocity, traj.angle_rad)
            traj.update(meteor, rho, dt)

            temperature_profile.append(meteor.core_temperature)
            t += dt

            # End if core sublimates (volatile, low pressure regime)
            if meteor.core_sublimated:
                sim_active = False
                result_text = "Core sublimated in thin atmosphere.\nSurvival chance: 0%"

        if traj.has_hit_ground() or meteor.is_fully_ablated() or t >= max_sim_time:
                sim_active = False

                if meteor.is_fully_ablated():
                    # Nothing survives if the meteor is completely destroyed
                    end_msg = "Meteor fully ablated in the atmosphere."
                    survival_percent = 0.0
                elif traj.has_hit_ground():
                    # Only compute survival if meteor reached ground
                    survival_prob = life.survival_probability(temperature_profile, dt)
                    survival_percent = round(100.0 * survival_prob, 1)
                    end_msg = "Meteor reached the ground."
                else:
                    end_msg = "Simulation time cap reached."
                    survival_prob = life.survival_probability(temperature_profile, dt)
                    survival_percent = round(100.0 * survival_prob, 1)

                result_text = f"{end_msg}\nSurvival chance for {life.name}: {survival_percent}%"

        # Render with time and layer
        ui.update_display(meteor, traj, atm, time_elapsed=t)
        if not sim_active and result_text is not None:
            ui.show_result_overlay(result_text)
        ui.flip()

    ui.shutdown()


if __name__ == "__main__":
    main()