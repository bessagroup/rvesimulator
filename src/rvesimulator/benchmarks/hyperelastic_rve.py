#                                                                       Modules
# =============================================================================
# Standard
import logging
import os
import time
from pathlib import Path
from typing import Any

import rvesimulator

from ..abaqus2py.abaqus_simulator import AbaqusSimulator
from ..microstructure.circle_particles import CircleParticles
from .py3rve_base import Py3RVEBase

#                                                         Authorship & Credits
# =============================================================================
__author__ = "Harikrishnan Vijayakumaran (H.Vijayakumaran@tudelft.nl)"
__credits__ = ["Harikrishnan Vijayakumaran", "Jiaxiang Yi"]
__status__ = "In progress"
# =============================================================================

# =============================================================================


class HyperelasticRVE(Py3RVEBase):
    """uni-axial tension for pp/pe composite without cohesive elements in
    between inclusion and matrix material phases"""

    def __init__(self) -> None:

        logging.basicConfig(level=logging.INFO,
                            filename='rve_simulation.log', filemode='w')
        self.logger = logging.getLogger("abaqus_simulation")
        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path":  Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "benchmark_abaqus_scripts.hyperelastic_rve",
            "sim_func": "HyperelasticRVE",
            "post_script": "benchmark_abaqus_scripts.hyperelastic_rve",
            "post_func": "PostProcess",
        }

    def update_sim_info(
        self,
        size: float = 0.02,
        radius_mu: float = 0.0031,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        params_matrix: list = [180.5, 2.8, 0.0025],
        params_inclusion: list = [950, 0.0005],
        mesh_division: int = 100,
        displacement_gradient: list = [[0.1, 0.0], [0.0, 0.0]],
        num_pseudo_time_steps: int = 10,
        simulation_time: float = 1.0,
        num_cpu: int = 8,
        model_name: str = "rve_2phase",
        part_name: str = "rve_2phase",
        instance_name: str = "rve_2phase",
        job_name: str = "hyperelastic_rve",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # micro_structure information
        self.seed = seed
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        # material properties
        self.params_matrix = params_matrix
        self.params_inclusion = params_inclusion
        # simulation information
        self.mesh_division = mesh_division
        self.displacement_gradient = displacement_gradient
        self.num_pseudo_time_steps = num_pseudo_time_steps
        self.simulation_time = simulation_time
        # model information
        self.model_name = model_name
        self.part_name = part_name
        self.instance_name = instance_name
        self.job_name = job_name
        # parallel information and platform
        self.num_cpu = num_cpu
        # get the micro_structure information
        self.seed = seed

        # update simulation information to logger
        self.logger.info("=============== simulation information ============")
        self.logger.info("size: {}".format(size))
        self.logger.info("radius_mu: {}".format(radius_mu))
        self.logger.info("radius_std: {}".format(radius_std))
        self.logger.info("vol_req: {}".format(vol_req))
        self.logger.info("params_matrix: {}".format(params_matrix))
        self.logger.info("params_inclusion: {}".format(params_inclusion))
        self.logger.info("mesh_division: {}".format(mesh_division))
        self.logger.info("displacement_gradient: {}".format(displacement_gradient))
        self.logger.info(
            "num_pseudo_time_steps: {}".format(num_pseudo_time_steps))
        self.logger.info("simulation_time: {}".format(simulation_time))
        self.logger.info("num_cpu: {}".format(num_cpu))

        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "params_matrix": params_matrix,
            "params_inclusion": params_inclusion,
            "mesh_division": mesh_division,
            "displacement_gradient": displacement_gradient,
            "num_pseudo_time_steps": num_pseudo_time_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu}

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        self.sim_info = {
            "inclusion_location_information":
            self.microstructure.microstructure_info[
                "inclusion_location_information"],
            "radius_mu": self.microstructure.microstructure_info["radius_mu"],
            "radius_std": self.microstructure.microstructure_info[
                "radius_std"],
            "length_start":
            self.microstructure.microstructure_info["length_start"],
            "length_end":
            self.microstructure.microstructure_info["length_end"],
            "width_start":
            self.microstructure.microstructure_info["width_start"],
            "width_end": self.microstructure.microstructure_info["width_end"],
            "model_name": self.model_name,
            "part_name": self.part_name,
            "instance_name": self.instance_name,
            "job_name": self.job_name,
            "displacement_gradient": self.displacement_gradient,
            "params_matrix": self.params_matrix,
            "params_inclusion": self.params_inclusion,
            "mesh_division": self.mesh_division,
            "num_pseudo_time_steps": self.num_pseudo_time_steps,
            "simulation_time": self.simulation_time,
            "num_cpu": self.num_cpu
            }

    def run_simulation(
        self,
        sample: dict = None,
        folder_index: int = None,
        delete_odb: bool = True,
    ) -> dict:
        """run single simulation

        Parameters
        ----------
        sample : dict, optional
            a dict contains the information of design variables
        folder_index : int, optional
            first folder index, by default None
        delete_odb : bool, optional

        Returns
        -------
        dict
            all the simulation results from abaqus
        """
        # number of samples
        self._create_working_folder(folder_index)

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
            simulator.run(py_func=self.folder_info["sim_func"],
                          py_script=self.folder_info["sim_script"],
                          post_py_func=self.folder_info["post_func"],
                          num_cpu=self.num_cpu,
                          post_py_script=self.folder_info["post_script"],
                          delete_odb=delete_odb)
            # get the simulation results back
            results = simulator.read_back_results()
            self.logger.info("abaqus simulation finished")
        except FileNotFoundError:
            self.logger.error("simulation failed")
            results = None
        end_time = time.time()
        self.logger.info("time used: {} s".format(end_time - start_time))
        self.logger.info("============== End abaqus simulation ============")

        # back to main folder
        os.chdir(self.main_folder)

        return results
