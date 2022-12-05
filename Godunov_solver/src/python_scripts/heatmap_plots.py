from dash import Dash, dcc, html, Input, Output
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly import tools
import pandas as pd
import sys
import numpy as np
import plotly.express as px

n_cells = int(sys.argv[1])
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

data = pd.read_csv(data_dir + 'piston_r_u_p.txt', delimiter='\t', names=['r', 'u', 'p', 't'])
data['x'] = data.index
length = int(len(data.index) / n_cells)

GAMMA = 5. / 3
c_v = 8.31 / (GAMMA - 1)
data['s'] = 1 + c_v * np.log(data['p'] / (data['r'] ** GAMMA))


def data_step(i):
    return data.iloc[i * n_cells: (i + 1) * n_cells]


times = []
for step in np.arange(0, length - 1):
    times.append(str(data_step(step)["t"].iloc[0]))

heat_data = []

# for i in range(0, length - 1):
#     for j in range(len(data_step(i).x)):
#         t = data_step(i)["t"].iloc[0]
#         x = data_step(i).x.iloc[j]
#         r = data_step(i).r.iloc[j]
#         print(t, x, r)
#         heat_data.append([t, x, r])
#
# fig = px.imshow(heat_data)
# fig.show()


dx = dy = 0.05
y, x = np.mgrid[-5 : 5 + dy : dy, -5 : 10 + dx : dx]
z = np.sin(x)**10 + np.cos(10 + y*x) + np.cos(x) + 0.2*y + 0.1*x

print(z)