
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
    seed = int(design['seed'])
               

    task = HyperelasticRVE()
    task.update_sim_info(size=4,
                         radius_mu=0.15,
                         radius_std=0.0,
                         params_matrix=[3.025e-1, 5.083, 2.617e-1],
                         params_inclusion=[2.645e02, 1.745e-3],
                         mesh_division=300,
                         platform='cluster',
                         vol_req=vol_req,
                         num_cpu=12,
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
