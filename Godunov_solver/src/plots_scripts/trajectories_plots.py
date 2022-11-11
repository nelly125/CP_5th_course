import pandas as pd
import matplotlib.pyplot as plt
import sys

dir_name = sys.argv[2]

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
# ax4 = fig.add_subplot(224)

for a in [ax1, ax2, ax3]:
    for label in (a.get_xticklabels() + a.get_yticklabels()):
        label.set_fontsize(16)

trajectories_data = pd.read_csv(dir_name + "trajectories.txt", delimiter="\t", names = ["t", "left", "diaph", "right"])

ax1.plot(trajectories_data["t"], trajectories_data["left"], color='red')
ax1.set_title("Left side", fontsize=20)
ax2.plot(trajectories_data["t"], trajectories_data["right"], color='green')
ax2.set_title("Right side", fontsize=20)
ax3.plot(trajectories_data["t"], trajectories_data["diaph"], color='purple')
ax3.set_title("Contact trajectory", fontsize=20)

# plt.tick_params(axis='both', which='major', labelsize=17)

fig.savefig(dir_name + "plots/" + "trajectories.png")
plt.show(block=False)
plt.close()