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
from rvesimulator.design_of_experiment.path_sampler import StrainPathSampler
from rvesimulator.design_of_experiment.samplers import FixNumberSampler

# design variable is set to be the number of control points for loads path
# therefore the
doe_variables = OrderedDict({"num_control": 7, "num_increment": 100})

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
# print(data)

# initialize the strain path generator
strain_path_generator = StrainPathSampler(seed=12, num_dim=6)
data = strain_path_generator.get_strain_path(
    data=data, arg_name="loads_path", interploation_method="quadratic"
)
strain_path_generator.plot_path(iteration=0)
