import sys
import os

dir_name = sys.argv[2]

n_cells = dir_name.split('_')[2]

# os.system("python3 ./src/python_scripts/animated_plots.py " + n_cells + " " + dir_name)

os.system("python3 ./src/python_scripts/plotly_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/contour_plots.py " + n_cells + " " + dir_name)
# os.system("python3 ./src/python_scripts/energy_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/delta_energy_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/trajectories_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/fourier_transform.py " + n_cells + " " + dir_name)
