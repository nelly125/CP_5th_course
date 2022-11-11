import sys

import matplotlib.pyplot as plt
import pandas as pd

N = sys.argv[1]
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

delta_energy_data = pd.read_csv(data_dir + 'delta_energy.txt', delimiter='\t',
                                names=['t', 'delta_energy', 'A_l', 'A_r'])

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)

for a in [ax1, ax2, ax3]:
    for label in (a.get_xticklabels() + a.get_yticklabels()):
        label.set_fontsize(16)

ax1.plot(delta_energy_data["t"], delta_energy_data["delta_energy"], color='red')
ax1.set_title("delta energy", fontsize=20)
ax2.plot(delta_energy_data["t"], delta_energy_data["A_l"], color='green')
ax2.set_title("A left piston", fontsize=20)
ax3.plot(delta_energy_data["t"], delta_energy_data["A_r"], color='purple')
ax3.set_title("A right piston", fontsize=20)

fig.savefig(dir_name + "plots/" + "delta_energy.png")
plt.show(block=False)
plt.close()
