# import needed libraries
from rvesimulator.benchmarks.hyperelastic_rve import HyperelasticRVE

import numpy as np
import pandas as pd


def run_simulation():
    # load dataframe from csv
    displacement_gradient = pd.read_csv('input_data_space.csv')

    # initialization
    task1 = HyperelasticRVE()
    # update simulation info
    task1.update_sim_info(size = 0.048,
                        radius_mu = 0.002,
                        radius_std = 0.0004,
                        vol_req= 0.4,
                        platform='cluster',
                        seed=17)

    task1_results = {}
    # calculate responses of simulation
    for ii in range(len(displacement_gradient)):
        task1_results[ii] = task1.run_simulation(
            sample={'displacement_gradient': displacement_gradient.iloc[ii].values.reshape(2,2).tolist()}, 
            third_folder_index=ii
        )

def main():
    run_simulation()