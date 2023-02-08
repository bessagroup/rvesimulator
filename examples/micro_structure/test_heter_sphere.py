# third-party packages
import sys

# from collections import OrderedDict

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)
# local packages

from rvesimulator.microstructures.heter_radius_sphere import HeterRadiusSphere

# define the geometry information of the RVE
length = 1.0
width = 1.0
height = 1.0
radius = 0.1
radius_std = 0.02
vol_req = 0.15

CircleInclusionGenerator = HeterRadiusSphere(
    length=length,
    width=width,
    height=height,
    radius_mu=radius,
    radius_std=radius_std,
    vol_req=vol_req,
)


CircleInclusionGenerator.generate_rve()
CircleInclusionGenerator.plot_rve()
print(CircleInclusionGenerator.fiber_positions)
# CircleInclusionGenerator.save_results(file_name="3d_rve.json")
