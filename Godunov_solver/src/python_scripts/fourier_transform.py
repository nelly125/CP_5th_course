import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
import scipy

import sys

n_cells = int(sys.argv[1])
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

# data = pd.read_csv(data_dir + 'cell_params.txt', delimiter='\t', names=['t', 'r', 'u', 'p'], index_col=False)
# length = int(len(data.index) / n_cells)
#
# contact_cell = data.copy()
#
# param = 'p'
#
# p_mean = contact_cell[param].mean()
# contact_cell[param] -= p_mean
#
# p = contact_cell[param].iloc[:].to_numpy()
# t = contact_cell['t'].iloc[:].to_numpy()
# tr = np.linspace(min(t), max(t), len(t))
# vr = scipy.signal.resample(p, len(t))
#
# L = len(tr)
# Ts = np.mean(np.diff(tr))
# Fs = 1 / Ts
# Fn = Fs / 2
# FTvr = fftpack.fft(vr) / L
# Fv = np.linspace(0, 1, L // 2 + 1) * Fn
#
# fig = plt.figure(figsize=(18, 10))
# ax = fig.add_subplot(111)
# plt.plot(Fv, np.abs(FTvr[0:len(Fv)] * 2))
# ax.set_xlabel('Frequency', fontsize=16)
# ax.set_ylabel('Amplitude', fontsize=16)
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.title("Fourier transform", fontsize=20)
# fig.savefig(dir_name + "plots/" + "fourier_transform.png", transparent=False)
# plt.show(block=False)
# plt.close()
data = pd.read_csv(data_dir + 'piston_r_u_p.txt', delimiter='\t', names=['r', 'u', 'p', 't'])
data['x'] = data.index

n_steps = int(len(data) / n_cells)

point = 2.5

p_i= []
for i in range(n_steps):
    x = data.iloc[i*n_cells:(i + 1)*n_cells]["x"]
    p = data.iloc[i*n_cells:(i + 1)*n_cells]["p"]
    p_new = np.interp(point, x, p)
    p_i.append(p_new)

t_i = data["t"].unique()
contact_cell = pd.DataFrame(p_i, columns=["p"])
contact_cell["t"] = t_i
contact_cell = contact_cell.iloc[:len(contact_cell) // 2]
p_mean = contact_cell['p'].mean()
contact_cell['p'] -= p_mean

p = contact_cell['p'].iloc[:].to_numpy()

t = contact_cell['t'].iloc[:].to_numpy()

tr = np.linspace(min(t), max(t), len(t))
vr = scipy.signal.resample(p, len(t))

L = len(tr)
Ts = np.mean(np.diff(tr))
Fs = 1/Ts
Fn = Fs/2
FTvr = fftpack.fft(vr)/L
Fv = np.linspace(0, 1, L//2 +1) * Fn
len(FTvr[1:len(Fv)])

fig = plt.figure(figsize=(18, 10))
ax = fig.add_subplot(111)
plt.plot(Fv, np.abs(FTvr[0:len(Fv)] * 2))
ax.set_xlabel('Frequency', fontsize=16)
ax.set_ylabel('Amplitude', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(0, 5)
plt.title("Fourier transform", fontsize=20)
fig.savefig(dir_name + "plots/" + "fourier_transform.png", transparent=False)
plt.show(block=False)
plt.close()