# third-party packages
import sys

# from collections import OrderedDict

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)
# local packages

from rvesimulator.microstructures.heter_radius_circles import (
    HeterCircleInclusion,
)

# define the geometry information of the RVE
length = 0.048
width = 0.048
radius = 0.01
radius_std = 0.005
vol_req = 0.45

CircleInclusionGenerator = HeterCircleInclusion(
    length=length,
    width=width,
    radius_mu=radius,
    radius_std=radius_std,
    vol_req=vol_req,
)


CircleInclusionGenerator.generate_rve()
CircleInclusionGenerator.plot_rve(save_figure=False, fig_name="rve.png")
print(CircleInclusionGenerator.fiber_positions)
CircleInclusionGenerator.save_results()
