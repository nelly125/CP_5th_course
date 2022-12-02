import os
import re

import matplotlib.pyplot as plt
import pandas as pd

# dir_to_output = "./results/plots_to_report/"

# files_dir = 'results/amplitude_dependency/time_50/'
# files_dir = 'results/omega_dependency/'
# files_dir = 'results/sigma_dependency/'
files_dir = 'results/l_dependency/'

dirs = [f for f in os.listdir(files_dir) if re.match(r'piston*', f)]
# print(dirs)
i = 0
for dir in dirs:
    path = files_dir + dir + "/data/"
    parameters = dir.split("__")
    data_i = pd.read_csv(path + "trajectories.txt", delimiter='\t', names=['r', 'left', 'contact', 'right']).tail(1)
    data_i["l"] = float(parameters[6])
    data_i=data_i.drop(["r"], axis=1)
    if (i == 0):
        data = data_i
    else:
        pass
        data = pd.concat([data, data_i])
    i += 1


data=data.reindex(columns=["l", "left", "contact", "right"])
data=data.sort_values(by="l", ascending=True)
print(data)

data.to_csv(files_dir + "trajectories_comparison.txt", index=False, sep="\t")