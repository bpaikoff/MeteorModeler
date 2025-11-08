# validate.py
import math
from equations import Equations
from meteor import Meteor
from trajectory import Trajectory
from atmosphere import Atmosphere

def run_case(composition="ice", diameter=1.0, v0=19000.0, angle_deg=45.0, dt=0.02):
    eq = Equations()
    atm = Atmosphere()
    meteor = Meteor(composition=composition, diameter=diameter, equations=eq)
    traj = Trajectory(initial_altitude=120_000.0, initial_velocity=v0, entry_angle_deg=angle_deg)

    t = 0.0
    results = []

    while not traj.has_hit_ground() and not meteor.is_fully_ablated() and not meteor.core_sublimated and t < 600.0:
        rho = atm.get_density(traj.altitude)
        T_air = atm.get_temperature(traj.altitude)
        P_air = atm.get_pressure(traj.altitude)
        Pvap = eq.vapor.vapor_pressure("ice", meteor.temperature)
        #print(f"t={t:.2f}s, T={meteor.temperature:.1f}K, Pvap={Pvap:.2e} Pa, Pamb={P_air:.2e} Pa")

        eq.update_thermo_ablation(meteor, rho, T_air, traj.velocity, dt, P_air, traj.angle_rad)
        eq.check_breakup(meteor, rho, traj.velocity, traj.angle_rad)
        traj.update(meteor, rho, dt)

        results.append({
            "time": t,
            "alt_km": traj.altitude / 1000.0,
            "vel": traj.velocity,
            "outer_T": meteor.temperature,
            "core_T": meteor.core_temperature,
            "phase": meteor.phase,
            "mass": meteor.mass,
        })
        t += dt

    return results

if __name__ == "__main__":
    # Example: run ice meteoroid case
    data = run_case(composition="ice", diameter=1.0, v0=19000.0, angle_deg=45.0)
    for row in data[::50]:  # print every 50th step
        print(f"{row['time']:6.2f}s | {row['alt_km']:6.1f} km | v={row['vel']:7.0f} m/s | "
              f"T_out={row['outer_T']:6.0f}K | T_core={row['core_T']:6.0f}K | "
              f"phase={row['phase']} | mass={row['mass']:8.1f} kg")