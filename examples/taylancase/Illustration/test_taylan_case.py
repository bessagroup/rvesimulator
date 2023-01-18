# import needed libraries
import os
import pickle
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)

# import local packages and functions
import rvesimulator
from rvesimulator.design_of_experiment.samplers import FixNumberSampler
from rvesimulator.simulators.path_generator import PathGenerator
from rvesimulator.simulators.taylan_rve import TaylanRVE

# design of experiments
# here the design variable is set to be the number of control points for loads path
# therefore the
doe_variables = OrderedDict({"control_points": 7})

# define number of samples, which means three different path will be generated
num_points = 3

# define the information of outputs
name_outputs = ["stress", "strain", "plastic_energy"]
doe_sampler = FixNumberSampler()
doe_sampler.sampling(
    num_samples=num_points,
    design_space=doe_variables,
    out_names=name_outputs,
    seed=123456,
)
data = doe_sampler.data
data = doe_sampler.data
print(doe_sampler.data)

# create an empty dataframe that contain the keys of loads_path
loads_path_temp = np.empty([num_points, 1])
loads_path_temp[:] = np.nan
loads_path = pd.DataFrame(loads_path_temp, columns=["loads_path"])
loads_path["loads_path"] = loads_path["loads_path"].astype(object)

# define the path generator
for ii in range(num_points):
    path_generator = PathGenerator(
        num_control_points=data["samples"].at[ii, "control_points"]
    )
    loads_path.iloc[ii, 0] = path_generator.quadratic_interpolate()
    # path_generator.plot_path()

# add the loads path to the samples
data["samples"] = pd.concat(
    [data["samples"], loads_path], axis=1, join="inner"
)
print(data["samples"])

# load the abaqus simulation modulus and run the simulation
simulation_wrapper = TaylanRVE()
simulation_wrapper.update_sim_info(
    print_info=True, mesh_partition=100, num_cpu=1, task="task1"
)

simulation_wrapper.run_simulation(data=data)
