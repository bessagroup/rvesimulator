# third-party packages
import sys

# from collections import OrderedDict

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)
# local packages

from matplotlib import pyplot as plt

from rvesimulator.microstructures.heter_radius_circles import (
    HeterCircleInclusion,
)

# define the geometry information of the RVE
length = 1
width = 1
radius = 0.1
radius_std = 0.02
vol_req = 0.40

CircleInclusionGenerator = HeterCircleInclusion(
    length=length,
    width=width,
    radius_mu=radius,
    radius_std=radius_std,
    vol_req=vol_req,
    seed=1,
)


CircleInclusionGenerator.generate_rve()
CircleInclusionGenerator.plot_rve(save_figure=False, fig_name="rve.png")
# print(CircleInclusionGenerator.fiber_positions)
CircleInclusionGenerator.save_results()
rgmsh = CircleInclusionGenerator.crate_rgmsh(num_discrete=600)
CircleInclusionGenerator.rgmsh_plot()
