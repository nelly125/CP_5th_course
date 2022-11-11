import sys
import os

n_cells = sys.argv[1]
dir_name = sys.argv[2]

# os.system("python3 ./src/python_scripts/animated_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/energy_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/delta_energy_plots.py " + n_cells + " " + dir_name)
os.system("python3 ./src/python_scripts/trajectories_plots.py " + n_cells + " " + dir_name)
