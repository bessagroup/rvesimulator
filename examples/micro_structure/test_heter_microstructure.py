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
length = 1.0
width = 1.0
radius = 0.04
radius_std = 0.004
vol_req = 0.40

CircleInclusionGenerator = HeterCircleInclusion(
    length=length,
    width=width,
    radius_mu=radius,
    radius_std=radius_std,
    vol_req=vol_req,
)

#
for ii in range(1):
    CircleInclusionGenerator.generate_rve()
    CircleInclusionGenerator.plot_rve(
        save_figure=True, fig_name="rve_" + str(ii) + ".png"
    )
print(CircleInclusionGenerator.fiber_positions[:, 3])
