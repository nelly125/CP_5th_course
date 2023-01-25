import sys
import os

dir_name = sys.argv[1]

n_cells = dir_name.split('_')[2]

# os.system("python3 ./src/python_scripts/animated_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/plotly_plots_spherical.py " + n_cells + " " + dir_name)