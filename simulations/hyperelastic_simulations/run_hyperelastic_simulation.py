# import needed libraries
from rvesimulator.benchmarks.hyperelastic_rve import HyperelasticRVE

import numpy as np
import pandas as pd


def run_simulation():
    # load dataframe from csv
    displacement_gradient = pd.read_csv('doe.csv')

    # initialization
    task = HyperelasticRVE()
    # update simulation info
    task.update_sim_info(size=4,
                         radius_mu=0.2,
                         radius_std=0.04,
                         vol_req=0.4,
                         mesh_division=500,
                         platform='cluster',
                         seed=1234)

    task_results = pd.DataFrame(columns=['F11', 'F12', 'F21', 'F22', 'P11', 'P12', 'P21', 'P22'])

    # calculate responses of simulation
    for ii in range(len(displacement_gradient)):

        results = task.run_simulation(
            sample={'displacement_gradient': displacement_gradient.iloc[ii].values.reshape(2, 2).tolist()},
            third_folder_index=ii)
        required_result_keys = ['deformation_gradient', 'pk1_stress']
        if set(required_result_keys).issubset(results.keys()):
            required_results = np.concatenate([results[key][-1].reshape(-1, 4) for key in required_result_keys], axis=1).reshape(-1)
            required_keys = ['F11', 'F12', 'F21', 'F22', 'P11', 'P12', 'P21', 'P22']
            required_dict = dict(zip(required_keys, required_results))

            task_results = pd.concat(
                [task_results, pd.DataFrame([required_dict])], ignore_index=True)
            task_results.to_csv('output_data.csv', index=False)

if __name__ == "__main__":
    run_simulation()
