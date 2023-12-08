import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap

def create_discrete_cmap(name, color_list, bounds):
    cmap = mpl.colors.LinearSegmentedColormap.from_list(name, color_list, 256)
    cmap.set_under('steelblue')
    bounds = bounds
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return (cmap, norm)



yield_bounds = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0,
                85.0, 90.0, 95.0, 100.0]
hit_bounds = [0.0, 1.0]

yield_color_list =['#fffff5', '#edf7ea', '#dcefe0', '#cbe7d6', '#bbdfcd', '#acd6c4', '#9dcdbc', '#90c4b4', '#82bbad', '#76b2a6', '#6aa9a0', '#5e9f99', '#539694', '#488d8e', '#3e8389', '#337a84', '#29717f', '#1e677b', '#125e77', '#005573']

hit_color_list = ['#fcf9e2', '#000e50']

hit_cmap = create_discrete_cmap('hit_cmap', hit_color_list,hit_bounds)

yield_cmap = create_discrete_cmap('yield_cmap', yield_color_list, yield_bounds)
pass