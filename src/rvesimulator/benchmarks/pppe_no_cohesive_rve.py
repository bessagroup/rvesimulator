"""Class for uni-axial tension for pp/pe Mixture/composite without cohesive"""
#                                                                       Modules
# =============================================================================
# Standard
import logging
import os
import time
from typing import Any

# local
import rvesimulator

from ..abaqus2py.abaqus_simulator import AbaqusSimulator
from ..microstructure.circle_particles import CircleParticles
from .shared_functionalities import SimulationBase

#                                                         Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"

# =============================================================================

# =============================================================================


class PPPEMixtureNoCohesive(SimulationBase):
    """uni-axial tension for pp/pe Mixture/composite without cohesive elements
    in between fiber and matrix material phases"""

    def __init__(self) -> None:

        logging.basicConfig(level=logging.INFO,
                            filename='rve_simulation.log', filemode='w')
        self.logger = logging.getLogger("abaqus_simulation")
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(self.main_folder, "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                "/scriptbase/",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.pppe_mixture_no_coh",
            "sim_script": "PPPEMixtureNoCohesive",
            "post_path": "benchmark_abaqus_scripts.pppe_mixture_no_coh",
            "post_script": "PostProcess",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "benchmark_abaqus_scripts/vevp_leonov_model.f"

    def update_sim_info(
        self,
        size: float = 0.02,
        radius_mu: float = 0.0031,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            748.234,
            0.45,
            296.0,
            0.1,
            1.53e-02,
            6.608e-01,
            7.33e02,
            2.51e08,
            3.30e-26,
            7.42e06,
            2.67e-02,
            4.826e1,
            1.520,
            1.303,
        ],
        paras_pe: list = [
            611.6,
            0.45,
            296.0,
            0.1,
            4.40e-01,
            5.1e-03,
            772.1,
            2.29e08,
            9.94e-26,
            5.01e06,
            3.234e-01,
            15.88,
            13.52,
            0.043,
        ],
        mesh_partition: int = 100,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 1000,
        simulation_time: float = 100.0,
        num_cpu: int = 8,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        """update simulation information

        Parameters
        ----------
        size : float, optional
            size of rve, by default 0.02
        radius_mu : float, optional
            radius mean, by default 0.0031
        radius_std : float, optional
            radius standard deviation, by default 0.0
        vol_req : float, optional
            volume fraction requirement, by default 0.30
        paras_pp : list, optional
            parameters of pp
        paras_pe : list, optional
            parameters of pe
        mesh_partition : int, optional
            mesh partition, by default 100
        strain : list, optional
            maximum strain, by default [0.1, 0.0, 0.0]
        num_steps : int, optional
            number of simulation steps, by default 1000
        simulation_time : float, optional
            simulation time, by default 100.0
        num_cpu : int, optional
            number of cpu, by default 8
        platform : str, optional
            platform, by default "ubuntu"
        seed : Any, optional
            seed, by default None
        print_info : bool, optional
            print simulation information to the screen or not, by default False
        """
        # micro_structure information
        self.seed = seed
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        # material properties
        self.paras_pp = paras_pp
        self.paras_pe = paras_pe
        # simulation information
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        # parallel information and platform
        self.num_cpu = num_cpu
        self.platform = platform
        # get the micro_structure information
        self.seed = seed

        # update simulation information to logger
        self.logger.info("=============== simulation information ============")
        self.logger.info("size: {}".format(size))
        self.logger.info("radius_mu: {}".format(radius_mu))
        self.logger.info("radius_std: {}".format(radius_std))
        self.logger.info("vol_req: {}".format(vol_req))
        self.logger.info("paras_pp: {}".format(paras_pp))
        self.logger.info("paras_pe: {}".format(paras_pe))
        self.logger.info("mesh_partition: {}".format(mesh_partition))
        self.logger.info("strain: {}".format(strain))
        self.logger.info("num_steps: {}".format(num_steps))
        self.logger.info("simulation_time: {}".format(simulation_time))
        self.logger.info("num_cpu: {}".format(num_cpu))
        self.logger.info("platform: {}".format(platform))

        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform, }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        self.sim_info = {
            "job_name": "veni_nocoh_rve",
            "location_information": self.microstructure.microstructure_info[
                "location_information"
            ],
            "radius_mu": self.microstructure.microstructure_info["radius_mu"],
            "radius_std": self.microstructure.microstructure_info[
                "radius_std"],
            "len_start": self.microstructure.microstructure_info["len_start"],
            "len_end": self.microstructure.microstructure_info["len_end"],
            "wid_start": self.microstructure.microstructure_info["wid_start"],
            "wid_end": self.microstructure.microstructure_info["wid_end"],
            "paras_pp": self.paras_pp,
            "paras_pe": self.paras_pe,
            "mesh_partition": self.mesh_partition,
            "num_steps": self.num_steps,
            "simulation_time": self.simulation_time,
            "strain": self.strain,
            "num_cpu": self.num_cpu,
            "platform": self.platform,
            "subroutine_path": self.subroutine_path
        }

    def run_simulation(
        self,
        sample: dict = None,
        folder_index: int = None,
        sub_folder_index: int = None,
        third_folder_index: int = None,
    ) -> dict:
        """run single simulation

        Parameters
        ----------
        sample : dict, optional
            a dict contains the information of design variables
        folder_index : int, optional
            first folder index, by default None
        sub_folder_index : int, optional
            second folder index, by default None
        third_folder_index : int, optional
            third folder index, by default None

        Returns
        -------
        dict
            all the simulation results from abaqus
        """
        # number of samples
        self._create_working_folder(
            folder_index,
            sub_folder_index,
            third_folder_index,
        )
        os.chdir(self.working_folder)
        self.logger.info("working folder: {}".format(self.working_folder))
        # create microstructure
        self.microstructure = CircleParticles(
            length=self.size,
            width=self.size,
            radius_mu=self.radius_mu,
            radius_std=self.radius_std,
            vol_req=self.vol_req,
        )
        self.microstructure.generate_microstructure(seed=self.seed)
        self.microstructure.to_abaqus_format()
        self.microstructure.plot_microstructure(save_figure=True,
                                                fig_name="rve_{}.png".
                                                format(self.seed))
        self.vol_frac = self.microstructure.vol_frac
        self.logger.info("volume fraction: {}".format(self.vol_frac))
        # update simulation information
        self._get_sim_info()
        # update the geometry info for microstructure
        self._update_sample_info(sample=sample)
        # change folder to main folder
        # save microstructure
        # update simulation information
        self.logger.info("============== Start abaqus simulation ============")
        start_time = time.time()
        simulator = AbaqusSimulator(
            sim_info=self.sim_info, folder_info=self.folder_info
        )
        # run abaqus simulation
        try:
            simulator.execute()
            simulator.post_process(delete_odb=True)
            # get the simulation results back
            results = simulator.read_back_results()
            self.logger.info("simulation finished")
        except FileNotFoundError:
            # if the results file is not found, then we assume that there is an
            # error occurs in the simulation, return None and continue the
            # simulation
            self.logger.info("simulation failed")
            results = None
        end_time = time.time()
        self.logger.info("time used: {} s".format(end_time - start_time))
        self.logger.info("============== End abaqus simulation ============")

        # back to main folder
        os.chdir(self.main_folder)

        return results
