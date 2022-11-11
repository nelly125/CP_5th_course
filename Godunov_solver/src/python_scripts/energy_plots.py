import sys

import matplotlib.pyplot as plt
import pandas as pd

N = int(sys.argv[1])
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

energy_data = pd.read_csv(data_dir + 'energy.txt', delimiter='\t', names=['t', 'kinetic_energy', 'internal_energy'])
energy_data["total_energy"] = energy_data["kinetic_energy"] + energy_data["internal_energy"]

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)

for a in [ax1, ax2, ax3]:
    for label in (a.get_xticklabels() + a.get_yticklabels()):
        label.set_fontsize(16)

energy_data["total_energy"] = energy_data["total_energy"] / N
energy_data["kinetic_energy"] = energy_data["kinetic_energy"] / N
energy_data["internal_energy"] = energy_data["internal_energy"] / N

ax1.plot(energy_data["t"], energy_data["total_energy"], color='red')
ax1.set_title("Total energy", fontsize=20)
ax2.plot(energy_data["t"], energy_data["kinetic_energy"], color='green')
ax2.set_title("Kinetic energy", fontsize=20)
ax3.plot(energy_data["t"], energy_data["internal_energy"], color='purple')
ax3.set_title("Internal energy", fontsize=20)

fig.savefig(dir_name + "plots/" + "energy.png")
plt.show(block=False)
plt.close()
