import sys

import matplotlib.pyplot as plt
import pandas as pd

# data_dir_1 = "/Godunov_solver/results/amplitude_dependency/time_50/piston__1000__50.00__0.10__50.00__25.00__2.50/data/"
# data_dir_2 = "/Godunov_solver/results/amplitude_dependency/time_50/piston__1000__50.00__0.20__50.00__25.00__2.50/data/"
# data_dir_1 = "/home/nelly/Documents/MSU/git/CP_5th_course/Godunov_solver/results/omega_dependency/piston__1000__50.00__0.20__10.00__25.00__2.50/data/"
# data_dir_2 = "/home/nelly/Documents/MSU/git/CP_5th_course/Godunov_solver/results/omega_dependency/piston__1000__50.00__0.20__70.00__25.00__2.50/data/"

data_dir_1 = "/home/nelly/Documents/MSU/git/CP_5th_course/Godunov_solver/results/sigma_dependency/piston__1000__50.00__0.20__50.00__5.00__2.50/data/"
data_dir_2 = "/home/nelly/Documents/MSU/git/CP_5th_course/Godunov_solver/results/sigma_dependency/piston__1000__50.00__0.20__50.00__35.00__2.50/data/"


fig = plt.figure(figsize=(18, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)

for ax in [ax1, ax2, ax3]:
    ax.set_xlabel("s", fontsize=16)
    ax.set_ylabel("x", fontsize=16)
# ax4 = fig.add_subplot(224)

fig.tight_layout(pad=0.4, h_pad=4, w_pad=4)
plt.subplots_adjust(top=0.85, left=0.1, bottom=0.06)

for a in [ax1, ax2, ax3]:
    for label in (a.get_xticklabels() + a.get_yticklabels()):
        label.set_fontsize(16)

file = "trajectories.txt"

trajectories_data_1 = pd.read_csv(data_dir_1 + file, delimiter="\t", names=["t", "left", "diaph", "right"])
trajectories_data_2 = pd.read_csv(data_dir_2 + file, delimiter="\t", names=["t", "left", "diaph", "right"])

# first_name = "amplitude 0.1"
# second_name="amplitude 0.2"

first_name = "sigma 5"
second_name="sigma 35"

ax1.plot(trajectories_data_1["t"], trajectories_data_1["left"], color='red', label=first_name)
ax1.plot(trajectories_data_2["t"], trajectories_data_2["left"], color='blue', label=second_name)
ax1.set_title("Left side", fontsize=20)
ax2.plot(trajectories_data_1["t"], trajectories_data_1["right"], color='green', label=first_name)
ax2.plot(trajectories_data_2["t"], trajectories_data_2["right"], color='brown', label=second_name)
ax2.set_title("Right side", fontsize=20)
ax3.plot(trajectories_data_1["t"], trajectories_data_1["diaph"], color='purple', label=first_name)
ax3.plot(trajectories_data_2["t"], trajectories_data_2["diaph"], color='gold', label=second_name)
ax3.set_title("Contact trajectory", fontsize=20)

# plt.tick_params(axis='both', which='major', labelsize=17)

# fig.savefig(dir_name + "plots/" + "trajectories.png", transparent=False)
ax1.legend()
ax2.legend()
ax3.legend()
plt.show()
plt.close()


fig = plt.figure(figsize=(18, 10))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

for ax in [ax1, ax2]:
    ax.set_xlabel("s", fontsize=16)
    ax.set_ylabel("x", fontsize=16)
# ax4 = fig.add_subplot(224)

fig.tight_layout(pad=0.4, h_pad=4, w_pad=4)
plt.subplots_adjust(top=0.85, left=0.1, bottom=0.06)

for a in [ax1, ax2]:
    for label in (a.get_xticklabels() + a.get_yticklabels()):
        label.set_fontsize(16)

file = "energy.txt"

trajectories_data_1 = pd.read_csv(data_dir_1 + file, delimiter="\t", names=['t', 'kinetic_energy', 'internal_energy'])
trajectories_data_2 = pd.read_csv(data_dir_2 + file, delimiter="\t", names=['t', 'kinetic_energy', 'internal_energy'])

# first_name = "amplitude 0.1"
# second_name="amplitude 0.2"

ax1.plot(trajectories_data_1["t"], trajectories_data_1["kinetic_energy"], color='red', label=first_name)
ax1.plot(trajectories_data_2["t"], trajectories_data_2["kinetic_energy"], color='blue', label=second_name)
ax1.set_title("Kinetic_energy", fontsize=20)
ax2.plot(trajectories_data_1["t"], trajectories_data_1["internal_energy"], color='green', label=first_name)
ax2.plot(trajectories_data_2["t"], trajectories_data_2["internal_energy"], color='brown', label=second_name)
ax2.set_title("Internal", fontsize=20)

# plt.tick_params(axis='both', which='major', labelsize=17)

# fig.savefig(dir_name + "plots/" + "trajectories.png", transparent=False)
ax1.legend()
ax2.legend()
ax3.legend()
plt.show()
plt.close()