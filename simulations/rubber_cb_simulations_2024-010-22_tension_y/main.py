"""Tune hyper-parameters for each B.C. separately.
1. Load the dataset (of 25 problem BCs) - experimentdata for f3dasm
2. For each row, run optuna to tune the hyper-parameters
3. Collect all results and save to disk
"""
import logging
import os
import pickle
from pathlib import Path
from time import sleep

import f3dasm
import hydra
from f3dasm.logger import logger as f3dasm_logger
from run_hyperelastic_simulation import run_simulation  # run on each row!

f3dasm_logger.setLevel(logging.INFO)


def pre_process(config):
    """Loads the problem BC specifications as a dataset.
    This is used by f3dasm to run parallel jobs.

    Parameters
    ----------
    config
        Configuration parameters defined in config.yaml
    """
    if config.dataset.rerun:
        pickle_file_name = f"exp_data_{config.dataset.model}_jobs.pkl"
        job_status = pickle.load(open(Path.cwd() / pickle_file_name, 'rb'))
        job_status[job_status == 'in progress'] = 'open'
        job_status[job_status == 'error'] = 'open'
        # save the pickle again
        pickle.dump(job_status, open(Path.cwd() / pickle_file_name, 'wb'))
        return None
    # Actual dataset path
    dataset_path = Path(config.dataset.store_path_access)
    experimentdata = f3dasm.ExperimentData.from_csv(Path(
                            "{}/doe_rubber_cb_all_vf_tension_y.csv".format(dataset_path)))
    # Save to disk
    experimentdata.jobs.mark_all_open()
    experimentdata.store(filename='exp_data_{}'.format(config.dataset.model))
    return None


def process(config):
    """Process the data row-by-row to get the objective value.
    Each row is processed in parallel.

    Parameters
    ----------
    config
        Configuration parameters defined in config.yaml

    Returns
    -------
    data
        ExperimentData object. This is used by post_process to get the mean
        objective value.
    """
    data = f3dasm.ExperimentData.from_file(filename='exp_data_{}'.format(
                    'rve')) # Same name and location as before!
    print(os.getcwd())
    # call objective_per_example
    data.run(run_simulation, mode='cluster', # Change this to parallel if you want to run in parallel
             kwargs={'config': config}) # Run your function like this with any kwargs that you aht to pass on
    return None


@hydra.main(config_path="./", config_name="config", version_base="1.1")
def main_func(config):
    """ Main script distinguishes between the master and the workers."""
    # Optuna master - Runs and maintains the Job
    if config.hpc.jobid == 0: # Only the master executes this
        pre_process(config)
        process(config)
    # Workers
    else:
        print(config.hpc.jobid)
        sleep(3*config.hpc.jobid)  # Just to asynchronize the jobs
        process(config)


if __name__ == "__main__":
    main_func()
