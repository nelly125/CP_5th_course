import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

fig = plt.figure(figsize=(18, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)

trajectories_data = pd.read_csv(data_dir + "trajectories.txt", delimiter="\t", names=["t", "left", "diaph", "right"])

ax_values = {ax1: "left", ax2: "right", ax3: "diaph"}
ax_colors = {ax1: "red", ax2: "green", ax3: "purple"}
ax_names = {ax1: "Left boundary", ax2: "Right boundary", ax3: "Contact trajectory"}

fig.tight_layout(pad=6, h_pad=4, w_pad=6)
plt.subplots_adjust(top=0.9, left=0.1, bottom=0.06, wspace=0.1)

for ax in ax_values.keys():
    frate = 1. / 96
    ax.plot(trajectories_data["t"], trajectories_data[ax_values[ax]], color=ax_colors[ax])
    ax.set_title(ax_names[ax], fontsize=20)
    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('Amplitude', fontsize=16)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

fig.savefig(dir_name + "plots/" + "trajectories.png", transparent=False)
plt.show(block=False)
plt.close()
