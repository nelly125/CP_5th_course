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
c_v = 8.31 / (GAMMA - 1)
data['s'] = 1 + c_v * np.log(data['p'] / (data['r'] ** GAMMA))
data['Mach'] = data['u'] / np.sqrt(GAMMA * data['p'] / data['r'])


def data_step(i):
    return data.iloc[i * n_cells: (i + 1) * n_cells]


app = Dash(__name__)

fig = make_subplots(rows=2, cols=2)
for step in np.arange(0, length - 1):
    # fig.add_trace(go.Scatter(x=data_step(step).x, y=data_step(step).r,  name='density'), 1, 1)
    fig.add_trace(row=1, col=1, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="blue", width=3),
        name="density",
        x=data_step(step).x,
        y=data_step(step).r,

    ))
    fig.add_trace(row=1, col=2, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="green", width=3),
        name="velocity",
        x=data_step(step).x,
        y=data_step(step).u,
    ))
    fig.add_trace(row=2, col=1, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="red", width=3),
        name="pressure",
        x=data_step(step).x,
        y=data_step(step).p,
    ))
    fig.add_trace(row=2, col=2, trace=
    go.Scatter(
        visible=step == 0,
        line=dict(color="purple", width=3),
        name="Mach",
        x=data_step(step).x,
        y=data_step(step).Mach,
    ))
    # times.append(str(data_step(step)["t"].iloc[0]))

    # fig.add_trace(go.Scatter(x=data_step(step).x, y=data_step(step).u,  name='k(x)=cos(x)'), 1, 2)

window_length = 5

steps = []
for i in range(0, len(fig.data), 4):
    step = dict(
        method="restyle",
        # args=["visible", [False] * len(fig.data)],
        args=[
            # Make the ith trace visible
            {'visible': [False] * len(fig.data)},

            # Set the title for the ith trace
            {'title.text': 'Step %d' % i}],
        label=str(data_step(i // 4)["t"].iloc[0])
    )
    # print(step["args"][0]["visible"][i:i + 1])
    step["args"][0]["visible"][i:i + 1] = [True, True, True, True]
    steps.append(step)
    # fig.frames[]['layout'].update(title_text=f'My title {i}')

    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Time: ",
                      "suffix": " seconds",
                      "font": {"size": 25}},
        pad={"t": 50},
        steps=steps
    )]

title_font_size = 22
ticks_font_size = 20

fig['layout']['xaxis'].update(range=[data.x.min() - 0.05, data.x.max() + 0.05], title_text='x',
                              title_font={"size": title_font_size},
                              tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))
fig['layout']['yaxis'].update(title_text='density',
                              title_font={"size": title_font_size},
                              tickfont=dict(family='Rockwell', color='black', size=ticks_font_size),
                              type="log")
fig['layout']['xaxis2'].update(range=[data.x.min() - 0.05, data.x.max() + 0.05], title_text='x',
                               title_font={"size": title_font_size},
                               tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))
fig['layout']['yaxis2'].update(range=[data.u.min() - 0.05, data.u.max() + 0.05], title_text='velocity',
                               title_font={"size": title_font_size},
                               tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))
fig['layout']['xaxis3'].update(range=[data.x.min() - 0.05, data.x.max() + 0.05], title_text='x',
                               title_font={"size": title_font_size},
                               tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))
fig['layout']['yaxis3'].update(title_text='pressure',
                               title_font={"size": title_font_size},
                               tickfont=dict(family='Rockwell', color='black', size=ticks_font_size),
                               type="log")
fig['layout']['xaxis4'].update(range=[data.x.min() - 0.05, data.x.max() + 0.05], title_text='x',
                               title_font={"size": title_font_size},
                               tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))
fig['layout']['yaxis4'].update(range=[data.Mach.min() - 0.05, data.Mach.max() + 0.05], title_text='Mach',
                               title_font={"size": title_font_size},
                               tickfont=dict(family='Rockwell', color='black', size=ticks_font_size))

fig.update_layout(
    sliders=sliders,
)

fig.show()

fig.write_html(dir_name + "/plots/" + "animated_plot.html")

# chart_studio.tools.set_credentials_file(username='Nelly_125', api_key='IKyJ5ipPVE30cYXcOkyz')
# py.plot(fig, filename = 'file', auto_open=True)
