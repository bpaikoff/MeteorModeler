# montecarlo.py
import random
import matplotlib.pyplot as plt
from validate import run_case

compositions = ["iron", "stone", "ice", "carbonaceous"]
sizes = [0.1, 0.25, 0.5, 1.0, 2.0]  # meters
N_trials = 50  # number of random angles per composition/size

results = []

for comp in compositions:
    for d in sizes:
        for _ in range(N_trials):
            angle = random.uniform(20.0, 70.0)  # degrees
            data = run_case(composition=comp, diameter=d, v0=19000.0, angle_deg=angle)
            if data:
                peak_core_T = max(row["core_T"] for row in data)
                final_alt = data[-1]["alt_km"]
                results.append({
                    "composition": comp,
                    "diameter": d,
                    "angle": angle,
                    "peak_core_T": peak_core_T,
                    "final_alt": final_alt
                })

# Example visualization: scatter plot
for comp in compositions:
    xs = [r["diameter"] for r in results if r["composition"] == comp]
    ys = [r["peak_core_T"] for r in results if r["composition"] == comp]
    plt.scatter(xs, ys, label=comp)

plt.xlabel("Diameter (m)")
plt.ylabel("Peak Core Temperature (K)")
plt.legend()
plt.show()