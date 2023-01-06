# third-party packages
import sys

# from collections import OrderedDict

# path of local project
folder_path = "/home/yaga/Documents/GitHub/rvesimulator"
sys.path.insert(0, folder_path)
# local packages

from rvesimulator.microstructures.homo_radius_disks import CircleInclusion

# define the geometry information of the RVE
length = 1.0
width = 1.0
Radius = 0.05
Vol_req = 0.30

CircleInclusionGenerator = CircleInclusion(
    length=length,
    width=width,
    radius=Radius,
    vol_req=Vol_req,
    second_heuristic=False,
)
for ii in range(2):
    run_time = CircleInclusionGenerator.generate_rve()
    CircleInclusionGenerator.plot_rve(
        save_figure=True, fig_name="RVE_" + str(ii) + ".png"
    )
    CircleInclusionGenerator.save_results(
        file_name="test_" + str(ii) + ".json"
    )
