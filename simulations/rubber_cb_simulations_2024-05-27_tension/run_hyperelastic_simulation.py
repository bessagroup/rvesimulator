
import logging
import pickle
from pathlib import Path

import numpy as np
from f3dasm import Design
from f3dasm.logger import logger as f3dasm_logger

from rvesimulator.benchmarks.hyperelastic_rve import HyperelasticRVE


# IMPORTANT: the first argument has to be design,
# rest are the kwargs you passed
def run_simulation(design: Design, config):
    """Run one row of the dataset. Each BC 's hyper-parameters are tuned
    separately."""
    f3dasm_logger.info(f"*** Running on simulation #{design.job_number} ***")
    # Do what you want here
    # design has information about a row of the DOE you provided
    # can be accessed as design['INPUT_COLUMN_NAME'] to get the input value
    # Save outputs you want to f3dasm

    # get the design paramters of this row
    displacement_gradient = [[design['dU11'], design['dU12']], [design['dU21'], design['dU22']]]
    vol_req = design['vol_req']
    seed = 0
    
    # increase the possions ratio of the inclusion

               
    # reducing poisson's ratio for matrix to 0.3 and changing inclusion to 10%
    task = HyperelasticRVE()
    task.update_sim_info(size=4,
                         radius_mu=0.15,
                         radius_std=0.0,
                         params_matrix=[3.02495e-1, 5.08265, 6.67865e-1],  # abaqus takes mu as input, not mu_0
                         params_inclusion=[8.63832e1, 5.34292e-3],  # 500 times stiffer than matrix
                         num_pseudo_time_steps=20,
                         mesh_division=300,
                         platform='cluster',
                         vol_req=vol_req,
                         num_cpu=4,
                         seed=seed)
    
    try:
        task.run_simulation(sample={'displacement_gradient': displacement_gradient},
                            folder_index=design.job_number)

        # save status
        design['progress'] = 'finished'
        f3dasm_logger.info(f"*** simulation  \
                            #{design.job_number} has been finished ***")
    except:
        # save status
        design['progress'] = 'failed'
        f3dasm_logger.info(f"*** simulation  \
                            #{design.job_number} has been failed ***")
    return design
