import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import Slider, Button
from math import sqrt
import plotly.graph_objects as go


n_cells = int(sys.argv[1])
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

data = pd.read_csv(data_dir + 'piston_r_u_p.txt', delimiter='\t', names=['r', 'u', 'p', 't'])

length = int(len(data.index) / n_cells)

GAMMA = 5. / 3
data['x'] = data.index
data['Mach'] = data.г / sqrt(GAMMA * data.p/ data.r)


def data_step(i):
    return data.iloc[i * n_cells: (i + 1) * n_cells]


fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax1.set_xlabel('x', fontsize=16)
ax1.set_ylabel('density', fontsize=16)
ax1.set_xlim([data.x.min() - 0.05, data.x.max() + 0.05])
ax1.set_ylim([data.r.min() - 0.05, data.r.max() + 0.05])
# ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
# ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.grid(which='major')

ax2.set_xlabel('x', fontsize=16)
ax2.set_ylabel('velocity', fontsize=16)
ax2.set_xlim([data.x.min() - 0.05, data.x.max() + 0.05])
ax2.set_ylim([data.u.min() - 0.05, data.u.max() + 0.05])
# ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
# ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax2.grid(which='major')

ax3.set_xlabel('x', fontsize=16)
ax3.set_ylabel('pressure', fontsize=16)
ax3.set_xlim([data.x.min() - 0.05, data.x.max() + 0.05])
ax3.set_ylim([data.p.min() - 0.05, data.p.max() + 0.05])
# ax3.set_ylim([0.5 - 0.05, 2 + 0.05])
# ax3.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax3.grid(which='major')

ax4.set_xlabel('.x', fontsize=16)
ax4.set_ylabel('Mach', fontsize=16)
ax4.set_xlim([data.x.min() - 0.05, data.x.max() + 0.05])
ax4.set_ylim([data.Mach.min() - 0.05, data.Mach.max() + 0.05])
ax4.grid(which='major')

for a in [ax1, ax2, ax3]:
    for label in (a.get_xticklabels() + a.get_yticklabels()):
        label.set_fontsize(13)

ax1.set_title("time: " + str(data_step(0)["t"].iloc[0]) + " s")

density_plot, = ax1.plot(data_step(0).r, lw=2, color='blue', label="godunov_" + str(n_cells))
velocity_plot, = ax2.plot(data_step(0).u, lw=2, color='green', label="godunov_" + str(n_cells))
pressure_plot, = ax3.plot(data_step(0).p, lw=2, color='red', )
Mach_plt, = ax4.plot(data_step(0).Mach, lw=2, color='y', )

plt.subplots_adjust(left=0.25, bottom=0.25)
# plt.ion()

axfreq = plt.axes([0.25, 0.1, 0.65, 0.03])
freq_slider = Slider(
    ax=axfreq,
    label='Step',
    valmin=0,
    valmax=length - 1,
    valstep=1,
    valinit=0,
    color="purple"
)


def update(val):
    step = freq_slider.val
    density_plot.set_data(data_step(step).x, data_step(step).r)
    velocity_plot.set_data(data_step(step).x, data_step(step).u)
    ax1.set_title("time: " + str(data_step(step)["t"].iloc[0]) + " s")
    # ax1.autoscale()
    pressure_plot.set_data(data_step(step).x, data_step(step).p)
    Mach_plt.set_data(data_step(step).x, data_step(step).Mach)
    # ax3.autoscale()
    ax1.relim()
    ax2.relim()
    ax3.relim()
    ax4.relim()
    fig.canvas.draw_idle()


freq_slider.on_changed(update)
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    freq_slider.reset()


button.on_clicked(reset)

# ax2.legend()
# ax3.legend()

plt.show()
# plt.savefig("plot.,ng")
# plt.show(block=False)
# plt.pause(1)


plt.close()
