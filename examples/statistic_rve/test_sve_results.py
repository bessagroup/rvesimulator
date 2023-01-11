# import needed libraries
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)

from simulator_caller import SimulatorCaller

# import local packages and functions
import rvesimulator
from rvesimulator.design_of_experiment.samplers import FixNumberSampler
from rvesimulator.microstructures.microstructure_plots import DrawRVE2D
from rvesimulator.simulators.path_generator import PathGenerator

# create the doe
# define the design space
doe_variables = OrderedDict({"control_points": 7})

# define number of samples
num_points = 1
# define the information of outputs
name_outputs = ["vol_frac", "PK2", "Green_strain", "ALLPD"]
doe_sampler = FixNumberSampler()
doe_sampler.sampling(
    num_samples=num_points,
    design_space=doe_variables,
    out_names=name_outputs,
    seed=123456,
)
data = doe_sampler.data
print(data)

#
loads_path_temp = np.empty([num_points, 1])
loads_path_temp[:] = np.nan
# print(loads_path_temp)
loads_path = pd.DataFrame(loads_path_temp, columns=["loads_path"])
loads_path["loads_path"] = loads_path["loads_path"].astype(object)
# print(loads_path)

# define the path generator
for ii in range(num_points):
    path_generator = PathGenerator(
        num_control_points=data["samples"].at[ii, "control_points"]
    )
    loads_path.iloc[ii, 0] = path_generator.quadratic_interpolate()
    path_generator.plot_path()

# add the loads path to the samples
data["samples"] = pd.concat(
    [data["samples"], loads_path], axis=1, join="inner"
)
print(data)

simulation_wrapper = SimulatorCaller()
simulation_wrapper.update_sim_info(
    loads=[0.02, 0.02, 0.02], mesh_partition=100, print_info=True
)
result_bencmark_1 = simulation_wrapper.run_simulation(data=data)
