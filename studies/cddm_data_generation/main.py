#                                                                       Modules
# =============================================================================
# standard library

from pathlib import Path
from time import sleep

# third party
import f3dasm
import pandas as pd
from f3dasm import (CategoricalParameter, DiscreteParameter, Domain,
                    ExperimentData, ExperimentSample)
from f3dasm.datageneration import DataGenerator

from rvesimulator.additions.ampitudesampler import AmplitudeGenerator
# local functionalities
from rvesimulator.additions.hardening_law import (LinearHardeningLaw,
                                                  RambergHardeningLaw,
                                                  SwiftHardeningLaw)
from rvesimulator.benchmarks.cddm_rve import CDDM_RVE

# define the hardening law
swift_law = SwiftHardeningLaw(a=0.5, b=0.5, yield_stress=0.5)
#                                                          Authorship & Credits
# =============================================================================
__author__ = 'Jiaxiang Yi'
__credits__ = ['Jiaxiang Yi']
__status__ = 'Development'
# =============================================================================

class CDDMDataGenerator(DataGenerator):
    """

    Parameters
    ----------
    DataGenerator : DataGenerator
        data generator class from f3dasm
    """
    def execute(data: ExperimentSample) -> ExperimentSample:


        """execute simulation for CDDM data generator

        Parameters
        ----------
        data : ExperimentSample
            experimental data sample from f3dasm

        Returns
        -------
        ExperimentSample
            experimetal data with updated output data
        """
        # get one row of the experimental sample
        sample = data.experiment_sample
        # print job number to the screen
        print(sample.job_number)
        # get input parameters
        task_ID = sample.input_data['task_ID']
        # seed for generating amplitude curves
        seed = sample.input_data['seed']

        # initialize the ASCA RVE
        simulator = CDDMDataGenerator().initialize_simulator(task_ID=task_ID)

        # obtain amplitude curves
        path_sampler = AmplitudeGenerator(num_dim=3)
        paths = path_sampler.get_amplitude(
            num_amplitude=1,
            num_control=8,
            num_steps=100,
            arg_name="strain_amplitude",
            seed=seed,
        )
        samples_dict = paths.to_dict("records")
        # add output parameters
        try:
            simulator.run_simulation(sample=samples_dict[0],
                                     folder_index=sample.job_number)
            # get the output data
            sample.output_data['progress'] = "finished"

        except Exception:
            # if the job failed, then we set the progress to be failed
            sample.output_data['progress'] = "failed"

        return data

    def initialize_simulator(self, task_ID: int) -> None:

        if task_ID == "task_A":
            return self._config_simulator_of_task_A()
        elif task_ID == "task_B":
            return self._config_simulator_of_task_B()
        elif task_ID == "task_C":
            return self._config_simulator_of_task_C()
        elif task_ID == "task_D":
            return self._config_simulator_of_task_D()
        else:
            raise ValueError("task ID is not defined")

    def _config_simulator_of_task_A(self) -> None:

        # initialize CDDM RVE simulator
        simulator = CDDM_RVE()
        # update simulation parameters
        simulator.update_sim_info(
            mesh_partition=100,
            strain=[0.02, 0.02, 0.02],
            vol_req=0.45,
            radius_mu=0.01,
            radius_std=0.003,
            youngs_modulus_fiber=10,
            youngs_modulus_matrix=100,
            hardening_law=LinearHardeningLaw(a=0.5, yield_stress=0.5),
            num_cpu=8,
            platform="cluster",
            seed=2,
            print_info=False)
        return simulator

    def _config_simulator_of_task_B(self) -> None:
        # initialize CDDM RVE simulator
        simulator = CDDM_RVE()
        # update simulation parameters
        simulator.update_sim_info(
            mesh_partition=100,
            strain=[0.02, 0.02, 0.02],
            vol_req=0.30,
            radius_mu=0.003,
            radius_std=0.0,
            youngs_modulus_fiber=1,
            youngs_modulus_matrix=100,
            hardening_law=SwiftHardeningLaw(a=0.5, b=0.4, yield_stress=0.5),
            num_cpu=8,
            platform="cluster",
            seed=17,
            print_info=False)

        return simulator

    def _config_simulator_of_task_C(self) -> None:
        # initialize CDDM RVE simulator
        simulator = CDDM_RVE()
        # update simulation parameters
        simulator.update_sim_info(
            mesh_partition=100,
            strain=[0.02, 0.02, 0.02],
            vol_req=0.15,
            radius_mu=0.0015,
            radius_std=0.0003,
            youngs_modulus_fiber=1000,
            youngs_modulus_matrix=100,
            hardening_law=RambergHardeningLaw(a=0.5, b=0.4, yield_stress=0.5),
            num_cpu=8,
            platform="cluster",
            seed=23,
            print_info=False)

        return simulator

    def _config_simulator_of_task_D(self) -> None:
        # initialize CDDM RVE simulator
        simulator = CDDM_RVE()
        # update simulation parameters
        simulator.update_sim_info(
            mesh_partition=100,
            strain=[0.02, 0.02, 0.02],
            vol_req=0.30,
            radius_mu=0.003,
            radius_std=0.0,
            youngs_modulus_fiber=10,
            youngs_modulus_matrix=100,
            hardening_law=SwiftHardeningLaw(a=0.5, b=0.4, yield_stress=0.5),
            num_cpu=8,
            platform="cluster",
            seed=17,
            print_info=False)

        return simulator



def create_experimentdata() -> None:
    """this function is used to create design of experiment for CDDM case, the
    DOE contains two columns, one is the task ID number, the other is the seed
    used for generating the amplitude curves."""

    # TODO: read from file can not work properly from f3dasm
    # get the working path
    file_path = Path.cwd()
    # joint the path
    file_path = file_path.joinpath('design_of_experiments.csv')
    # read csv file via pandas
    samples = pd.read_csv(file_path, index_col=0, header=0)
    # define a domain for the design of experiments
    domain = Domain()
    domain.add('task_ID',
               CategoricalParameter(['task_A', 'task_B', 'task_C', 'task_D']))
    domain.add('seed', DiscreteParameter(lower_bound=1, upper_bound=1000))

    # create experiment data
    data = ExperimentData(domain=domain)
    data.sample(sampler='random', n_samples=4000, seed=1)

    # fill input_data with correct values
    data.input_data.data['task_ID'] = samples['task_ID'].values
    data.input_data.data['seed'] = samples['seed'].values
    # add output parameters
    # data = f3dasm.ExperimentData.from_file(filename='exp_{}'.format('cddm'))
    data.add_output_parameter('progress')

    # save data
    data.store(filename='exp_{}'.format('cddm'))

def execute_experimentdata() -> None:

    # load data from file
    data = f3dasm.ExperimentData.from_file(filename='exp_{}'.format('cddm'))
    # run the function
    data.evaluate(CDDMDataGenerator(), mode='cluster')
    data.store(filename='exp_{}'.format('cddm'))


if __name__ == '__main__':
    # run the function
    """ Main script distinguishes between the master and the workers."""
    print(f3dasm.HPC_JOBID)
    if f3dasm.HPC_JOBID == 0:
        create_experimentdata()
        execute_experimentdata()
    elif f3dasm.HPC_JOBID > 0:
        sleep(f3dasm.HPC_JOBID)
        execute_experimentdata()
