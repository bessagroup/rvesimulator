"""
CDDM RVE case
"""
#                                                                       Modules
# =============================================================================
# Standard
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict

# local
import rvesimulator
from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
from rvesimulator.additions.hardening_law import (
    SwiftHardeningLaw, HardeningLaw)

from .py3rve_base import Py3RVEBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


class J2Plasticity_3D(Py3RVEBase):
    """Interface between python and abaqus of the 3D RVE case

    Parameters
    ----------
    SimulationBase : class
        base class for simulation
    """

    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""

        logging.basicConfig(level=logging.INFO, filename="Ti6Al4V.log")
        self.logger = logging.getLogger("abaqus_simulation")

        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "benchmark_abaqus_scripts.j2_plastic_one_phase_rve",
            "sim_func": "j2_plastic_3d_rve",
            "post_script": "benchmark_abaqus_scripts.j2_plastic_one_phase_rve",
            "post_func": "post_process",
        }

    def update_sim_info(
        self,
        size: float = 0.048,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_partition: int = 100,
        strain: list = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 4,
        hardening_law: HardeningLaw = SwiftHardeningLaw(),
        seed: Any = None,
        print_info: bool = False,
        regular_load: bool = False,
        strain_amplitude: list = None,
    ) -> None:
        """regular rve for Ti6Al4V 3D case

        Parameters
        ----------
        size : float, optional
            size of the rve, by default 0.048
        radius_mu : float, optional
            radius mean of the rve, by default 0.003
        radius_std : float, optional
            radius deviation of the rve, by default 0.0
        vol_req : float, optional
            volume fraction requirement, by default 0.30
        youngs_modulus_matrix : float, optional
            youngs modulus of the matrix material, by default 100.0
        poisson_ratio_matrix : float, optional
            poisson ratio of the matrix material, by default 0.3
        youngs_modulus_fiber : float, optional
            youngs modulus of the fiber material, by default 1.0
        poisson_ratio_fiber : float, optional
            poisson ratio of the fiber material, by default 0.19
        mesh_partition : int, optional
            mesh partition for the edges, by default 30
        strain : list, optional
            applied maximum strain, by default [0.1, 0.0, 0.0]
        num_steps : int, optional
            number of simulation steps, by default 100
        simulation_time : float, optional
            total simulation time, by default 1.0
        num_cpu : int, optional
            number of cpu used for simulation, by default 1
        hardening_law : Any, optional
            hardening law, by default SwiftHardeningLaw()
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get simulation information
        self.size = size
        self.radius_mu = 0.1*size
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.hardening_law = hardening_law
        self.seed = seed
        self.regular_load = regular_load
        self.strain_amplitude = strain_amplitude
        # get hardening law
        self.hardening_table = hardening_law.calculate_hardening_table()
        if self.regular_load:
            self.sim_paras = {
                "size": size,
                "radius_mu": self.radius_mu,
                "youngs_modulus": youngs_modulus,
                "poisson_ratio": poisson_ratio,
                "hardening_table": self.hardening_table,
                "mesh_partition": mesh_partition,
                "strain": strain,
                "num_steps": num_steps,
                "simulation_time": simulation_time,
                "num_cpu": num_cpu, }
        else:
            self.sim_paras = {
                "size": size,
                "radius_mu": self.radius_mu,
                "youngs_modulus": youngs_modulus,
                "poisson_ratio": poisson_ratio,
                "hardening_table": self.hardening_table,
                "mesh_partition": mesh_partition,
                "strain": strain,
                "num_steps": num_steps,
                "simulation_time": simulation_time,
                "num_cpu": num_cpu,
                "strain_amplitude": strain_amplitude}

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        if self.regular_load:
            self.sim_info = {
                "job_name": "J2_plasticity_3D",
                "radius_mu": self.radius_mu,
                "len_start": -self.radius_mu,
                "len_end": self.size+self.radius_mu,
                "wid_start": -self.radius_mu,
                "wid_end": self.size+self.radius_mu,
                "hei_start": -self.radius_mu,
                "hei_end": self.size+self.radius_mu,
                "youngs_modulus": self.youngs_modulus,
                "poisson_ratio": self.poisson_ratio,
                "mesh_partition": self.mesh_partition,
                "hardening_table": self.hardening_table,
                "num_steps": self.num_steps,
                "simulation_time": self.simulation_time,
                "strain": self.strain,
                "num_cpu": self.num_cpu,
            }
        else:
            self.sim_info = {
                "job_name": "J2_plasticity_3D",
                "radius_mu": self.radius_mu,
                "len_start": -self.radius_mu,
                "len_end": self.size+self.radius_mu,
                "wid_start": -self.radius_mu,
                "wid_end": self.size+self.radius_mu,
                "hei_start": -self.radius_mu,
                "hei_end": self.size+self.radius_mu,
                "youngs_modulus": self.youngs_modulus,
                "poisson_ratio": self.poisson_ratio,
                "mesh_partition": self.mesh_partition,
                "hardening_table": self.hardening_table,
                "num_steps": self.num_steps,
                "simulation_time": self.simulation_time,
                "strain": self.strain,
                "num_cpu": self.num_cpu,
                "strain_amplitude": self.strain_amplitude,
            }

    def run_simulation(
        self,
        sample: Dict = None,
        folder_index: int = None,
        delete_odb: bool = True,
    ) -> Dict:
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
        )
        self.logger.info("working folder: {}".format(self.working_folder))
        # update simulation information
        self._get_sim_info()
        # update the geometry info for microstructure
        self._update_sample_info(sample=sample)
        # update logger on samples
        self.logger.info("==============        update info      ============")
        self.logger.info("sample: {}".format(sample))
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
                          post_py_script=self.folder_info["post_script"],
                          num_cpu=self.num_cpu,
                          delete_odb=delete_odb)
            results = simulator.read_back_results()
            self.logger.info("abaqus simulation finished")
        except FileNotFoundError:
            self.logger.info("abaqus simulation failed")
            results = None
        # get the simulation results back
        end_time = time.time()
        self.logger.info("time used: {} s".format(end_time - start_time))
        self.logger.info("============== End abaqus simulation ============")

        # back to main folder
        os.chdir(self.main_folder)

        return results
