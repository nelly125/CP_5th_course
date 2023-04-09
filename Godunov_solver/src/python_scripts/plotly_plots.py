import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash
from plotly.subplots import make_subplots

n_cells = int(sys.argv[1])
dir_name = sys.argv[2]
data_dir = dir_name + "/data/"

data = pd.read_csv(data_dir + 'piston_r_u_p.txt', delimiter='\t', names=['r', 'u', 'p', 't'])
data['x'] = data.index
length = int(len(data.index) / n_cells)

GAMMA = 5. / 3
R = 8.314
c_v = 8.31 / (GAMMA - 1)
data['s'] = 1 + c_v * np.log(data['p'] / (data['r'] ** GAMMA))
data['Mach'] = data['u'] / np.sqrt(GAMMA * data['p'] / data['r'])
data['Temp'] = data['p'] / (data['r'] * R) / (GAMMA - 1)

print(data['Temp'])


def data_step(i):
    return data.iloc[i * n_cells: (i + 1) * n_cells]


data["log_r"] = np.log(data['r'])

app = Dash(__name__)

#  column_widths=[100, 100, 100, 100], row_heights=[500]

fig = make_subplots(rows=2, cols=3,
                    subplot_titles=("density", "velocity", "pressure", "Mach", "entropy", 'temperature'),
                    specs=[[{}, {}, {}],
                           [{}, {}, {}]],
                    horizontal_spacing=0.05)
for step in np.arange(0, length - 1):
    # fig.add_trace(go.Scatter(x=data_step(step).x, y=data_step(step).r,  name='density'), 1, 1)
    fig.add_trace(row=1, col=1, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="blue", width=3),
        name="density",
        x=data_step(step).x,
        y=data_step(step).r,
        showlegend=False,
    ))
    fig.add_trace(row=1, col=2, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="green", width=3),
        name="velocity",
        x=data_step(step).x,
        y=data_step(step).u,
        showlegend=False,
    ))
    fig.add_trace(row=1, col=3, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="red", width=3),
        name="pressure",
        x=data_step(step).x,
        y=data_step(step).p,
        showlegend=False,
    ))
    fig.add_trace(row=2, col=1, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="purple", width=3),
        name="Mach",
        x=data_step(step).x,
        y=data_step(step).Mach,
        showlegend=False,
    ))

    fig.add_trace(row=2, col=2, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="black", width=3),
        name="s",
        x=data_step(step).x,
        y=data_step(step).s,
        showlegend=False,
    ))

    fig.add_trace(row=2, col=3, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="orange", width=3),
        name="s",
        x=data_step(step).x,
        y=data_step(step).Temp,
        showlegend=False,
    ))
    # fig.add_trace(go.Scatter(x=data_step(step).x, y=data_step(step).u,  name='k(x)=cos(x)'), 1, 2)

window_length = 6

steps = []
for i in range(0, len(fig.data), 6):
    step = dict(
        method="restyle",
        # args=["visible", [False] * len(fig.data)],
        args=[
            # Make the ith trace visible
            {'visible': [False] * len(fig.data)},

            # Set the title for the ith trace
            {'title.text': 'Step %d' % i}],
        label=str(round(data_step(i // 6)["t"].iloc[0], 3))
    )
    # print(step["args"][0]["visible"][i:i + 1])
    step["args"][0]["visible"][i:i + 1] = [True, True, True, True, True, True]
    steps.append(step)
    # fig.frames[]['layout'].update(title_text=f'My title {i}')

    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Time: ",
                      "suffix": "",
                      "font": {"size": 25}},
        pad={"t": 50},
        steps=steps
    )]

title_font_size = 18
ticks_font_size = 16

ax_values = {"xaxis": "r", "xaxis2": "u", "xaxis3": "p", "xaxis4": "Mach", "xaxis5": "s", "xaxis6": "Temp"}
ax_y_values = {"xaxis": "yaxis", "xaxis2": "yaxis2", "xaxis3": "yaxis3", "xaxis4": "yaxis4", "xaxis5": "yaxis5",
               "xaxis6": "yaxis6"}
values_names = {"r": "log(density)", "u": "velocity", "p": "pressure", "s": "entropy", "Mach": "Mach",
                "Temp": "Temperature"}
for ax in ax_values:
    fig['layout'][ax].update(range=[data.x.min() - 0.05, data.x.max() + 0.05], title_text='x',
                             title_font={"size": title_font_size},
                             tickfont=dict(family='Rockwell', color='black',
                                           size=ticks_font_size), dtick=0.2
                             )
    fig['layout'][ax_y_values[ax]].update(range=[data[ax_values[ax]].min() - 0.05, data[ax_values[ax]].max() + 0.05],
                                          # title_text=values_names[ax_values[ax]],
                                          title_font={"size": title_font_size},
                                          tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))

fig.update_layout(
    sliders=sliders,
    autosize=False,
    width=1700,
    height=1200, )

fig.show()

fig.write_html(dir_name + "/plots/" + "animated_plot.html")

# chart_studio.tools.set_credentials_file(username='Nelly_125', api_key='IKyJ5ipPVE30cYXcOkyz')
# py.plot(fig, filename = 'file', auto_open=True)
