import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
import scipy

import sys

n_cells = int(sys.argv[1])
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

data = pd.read_csv(data_dir + 'cell_params.txt', delimiter='\t', names=['t', 'r', 'u', 'p'], index_col=False)
length = int(len(data.index) / n_cells)

contact_cell = data.copy()

param = 'p'

p_mean = contact_cell[param].mean()
contact_cell[param] -= p_mean

p = contact_cell[param].iloc[:].to_numpy()
t = contact_cell['t'].iloc[:].to_numpy()
tr = np.linspace(min(t), max(t), len(t))
vr = scipy.signal.resample(p, len(t))

L = len(tr)
Ts = np.mean(np.diff(tr))
Fs = 1 / Ts
Fn = Fs / 2
FTvr = fftpack.fft(vr) / L
Fv = np.linspace(0, 1, L // 2 + 1) * Fn

fig = plt.figure(figsize=(18, 10))
ax = fig.add_subplot(111)
plt.plot(Fv, np.abs(FTvr[0:len(Fv)] * 2))
ax.set_xlabel('Frequency', fontsize=16)
ax.set_ylabel('Amplitude', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title("Fourier transform", fontsize=20)
fig.savefig(dir_name + "plots/" + "fourier_transform.png", transparent=False)
plt.show(block=False)
plt.close()