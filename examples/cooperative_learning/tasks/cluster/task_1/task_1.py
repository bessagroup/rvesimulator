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
from rvesimulator.simulators.coopative_learning_case import CooperativeRVE


def main() -> None:
    # design of experiments
    # here the design variable is set to be the number of control points for loads path
    # therefore the
    doe_variables = OrderedDict({"num_control": 7, "num_increment": 100})

    # define number of samples
    num_points = 10
    # define the information of outputs
    name_outputs = ["strain", "stress", "plastic_energy"]
    doe_sampler = FixNumberSampler()
    doe_sampler.sampling(
        num_samples=num_points,
        design_space=doe_variables,
        out_names=name_outputs,
        seed=1,
    )
    data = doe_sampler.data

    # create an empty dataframe that contain the keys of loads_path
    loads_path_temp = np.empty([num_points, 1])
    loads_path_temp[:] = np.nan
    loads_path = pd.DataFrame(loads_path_temp, columns=["loads_path"])
    loads_path["loads_path"] = loads_path["loads_path"].astype(object)

    # define the path generator
    strain_path_generator = StrainPathSampler(seed=1, num_dim=3)
    data = strain_path_generator.get_strain_path(
        data=data.copy(),
        arg_name="loads_path",
        interploation_method="quadratic",
    )

    # load the abaqus simulation modulus and run the simulation
    simulation_wrapper = CooperativeRVE()
    simulation_wrapper.update_sim_info(
        mesh_partition=100,
        vol_req=0.45,
        radius_mu=0.01,
        radius_std=0.005,
        E_fiber=10,
        E_matrix=100,
        update_micro_structure=True,
        hardening_law="linear",
        yield_stress=0.5,
        a=0.5,
        num_cpu=6,
        seed=23,
        print_info=False,
    )
    data = simulation_wrapper.run_simulation(data=data)


if __name__ == "__main__":
    main()
