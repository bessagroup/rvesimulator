"""Hollow plate simulation case, with different material models
and path loading"""
#                                                                       Modules
# =============================================================================
# Standard
import logging
import os
import time
from pathlib import Path, PurePath
from typing import Any, Dict

# third-party
import numpy as np

# local
import rvesimulator
from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
from rvesimulator.additions.hardening_law import LinearHardeningLaw

from ..additions.utils import rve_microstructure_plot
from .py3rve_base import Py3RVEBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


# =============================================================================


class HollowPlateBase(Py3RVEBase):

    def _get_sim_info(self) -> None:
        """get simulation info
        """
        self.sim_info = self.sim_paras

    def run_simulation(
        self,
        sample: Dict = None,
        folder_index: int = 0,
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
        self._create_working_folder(folder_index)

        # log the working folder
        self.logger.info("working folder: {}".format(self.working_folder))

        # update simulation information
        self._get_sim_info()
        # update the geometry info for microstructure
        self._update_sample_info(sample=sample)
        self.logger.info("==============     updated samples    ============")
        self.logger.info("sample: {}".format(sample))
        # plot the microstructure
        rve_microstructure_plot(
            fibers=np.atleast_2d([self.sim_info["size"]/2,
                                  self.sim_info["size"]/2,
                                  self.sim_info["radius"]]),
            size=self.size,
            save_fig=True)
        # change folder to main folder
        # save microstructure
        # update simulation information
        self.logger.info("============== Start abaqus simulation ============")
        start_time = time.time()
        simulator = AbaqusSimulator(
            sim_info=self.sim_info, folder_info=self.folder_info
        )
        # run abaqus simulation (in one goal, without submit inp separately)
        simulator.run(py_func=self.folder_info["sim_func"],
                      py_script=self.folder_info["sim_script"],
                      post_py_func=self.folder_info["post_func"],
                      post_py_script=self.folder_info["post_script"],
                      num_cpu=self.num_cpu,
                      delete_odb=delete_odb)
        # get the simulation results back
        results = simulator.read_back_results()
        end_time = time.time()
        self.logger.info("time used: {} s".format(end_time - start_time))
        self.logger.info("============== End abaqus simulation ============")

        # back to main folder
        os.chdir(self.main_folder)

        return results


class ElasticRegularLoads(HollowPlateBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        logging.basicConfig(level=logging.INFO,
                            filename="hollow_plate_simulation.log",)
        self.logger = logging.getLogger("abaqus_simulation")
        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "benchmark_abaqus_scripts.hollow_plate_sve",
            "sim_func": "ElasticRegularLoads",
            "post_script": "basic_analysis_scripts.post_process",
            "post_func": "PostProcess2D",
        }

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        num_cpu: int = 1,
        print_info: bool = False,
    ) -> None:
        """update simulation information

        Parameters
        ----------
        size : float, optional
            size of rve , by default 1.0
        radius : float, optional
            radius of hollow plate , by default 0.2
        youngs_modulus : float, optional
            youngs modulus, by default 100.0
        poisson_ratio : float, optional
            poisson ratio, by default 0.3
        mesh_portion : int, optional
            mesh portion, by default 30
        strain : list, optional
            applied strain, by default [0.1, 0.0, 0.0]
        num_cpu : int, optional
            number of cpu for simulation, by default 1
        platform : str, optional
            platform of execution, by default "ubuntu"
        print_info : bool, optional
            print the simulation information, by default False
        """
        # get the simulation information
        self.size = size
        self.radius = radius
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_cpu = num_cpu

        # update simulation information to logger
        self.logger.info("============== Simulation information ============")
        self.logger.info("size: {}".format(self.size))
        self.logger.info("radius: {}".format(self.radius))
        self.logger.info("youngs_modulus: {}".format(self.youngs_modulus))
        self.logger.info("poisson_ratio: {}".format(self.poisson_ratio))
        self.logger.info("mesh_partition: {}".format(self.mesh_partition))
        self.logger.info("strain: {}".format(self.strain))
        self.logger.info("num_cpu: {}".format(self.num_cpu))

        # have the simulation paras
        self.sim_paras = {
            "job_name": "hollow_plate",
            "radius": radius,
            "size": size,
            "youngs_modulus": youngs_modulus,
            "poisson_ratio": poisson_ratio,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_cpu": num_cpu,
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

# =============================================================================


class VonMisesPlasticRegularLoads(HollowPlateBase):

    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        logging.basicConfig(level=logging.INFO,
                            filename="hollow_plate_plastic_simulation.log",)
        self.logger = logging.getLogger("abaqus_simulation")
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "benchmark_abaqus_scripts.hollow_plate_sve",
            "sim_func": "VonMisesPlasticRegularLoads",
            "post_script": "basic_analysis_scripts.post_process",
            "post_func": "PostProcess2D",
        }

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        hardening_law: Any = LinearHardeningLaw(),
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        print_info: bool = False,
    ) -> None:
        """update simulation information

        Parameters
        ----------
        size : float, optional
            size of rve, by default 1.0
        radius : float, optional
            radius of hollow plate, by default 0.2
        youngs_modulus : float, optional
            youngs modulus, by default 100.0
        poisson_ratio : float, optional
            poisson ratio, by default 0.3
        mesh_partition : int, optional
            mesh partition, by default 30
        strain : list, optional
            strain, by default [0.1, 0.0, 0.0]
        hardening_law : Any, optional
            hardening law, by default LinearHardeningLaw()
        num_steps : int, optional
            number of step, by default 100
        simulation_time : float, optional
            simulation time, by default 1.0
        num_cpu : int, optional
            number of cpu been used, by default 1
        platform : str, optional
            platform of simulation, by default "ubuntu"
        print_info : bool, optional
            print simulation info to screen, by default False
        """
        # get the simulation information
        self.size = size
        self.radius = radius
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.hardening_law = hardening_law
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu

        # get hardening table
        self.hardening_table = self.hardening_law.calculate_hardening_table()

        # update simulation information to logger
        self.logger.info("============== Simulation information ============")
        self.logger.info("size: {}".format(self.size))
        self.logger.info("radius: {}".format(radius))
        self.logger.info("young_modulus: {}".format(youngs_modulus))
        self.logger.info("poisson_ratio: {}".format(poisson_ratio))
        self.logger.info("mesh_partition: {}".format(mesh_partition))
        self.logger.info("strain: {}".format(strain))
        self.logger.info("hardening_law: {}".format(hardening_law))
        self.logger.info("hardening_table: {}".format(self.hardening_table))
        self.logger.info("num_steps: {}".format(num_steps))
        self.logger.info("simulation_time: {}".format(simulation_time))
        self.logger.info("num_cpu: {}".format(num_cpu))

        # have the simulation paras

        self.sim_paras = {"job_name": "hollow_plate",
                          "radius": radius,
                          "size": size,
                          "youngs_modulus": youngs_modulus,
                          "poisson_ratio": poisson_ratio,
                          "hardening_table": self.hardening_table,
                          "mesh_partition": mesh_partition,
                          "strain": strain,
                          "num_steps": num_steps,
                          "simulation_time": simulation_time,
                          "num_cpu": num_cpu, }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)


# =============================================================================
class VonMisesPlasticPathLoads(HollowPlateBase):
    """Hollow plate simulation case, with von mises plasticity and path loading

    Parameters
    ----------
    HollowplateBase : class
        base class of hollow plate simulation case
    """

    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        # define the logger
        logging.basicConfig(level=logging.INFO,
                            filename="plastic_path_loads_simulation.log",)
        self.logger = logging.getLogger("abaqus_simulation")

        # folder information
        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir":  Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "benchmark_abaqus_scripts.hollow_plate_sve",
            "sim_func": "VonMisesPlasticPathLoads",
            "post_script": "basic_analysis_scripts.post_process",
            "post_func": "PostProcess2D",
        }

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        strain_amplitude: list = None,
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        hardening_law: Any = LinearHardeningLaw(),
        print_info: bool = False,
    ) -> None:
        """update simulation information

        Parameters
        ----------
        size : float, optional
            size of RVE, by default 1.0
        radius : float, optional
            radius of hollow plate, by default 0.2
        youngs_modulus : float, optional
            young modulus, by default 100.0
        poisson_ratio : float, optional
            poisson ratio, by default 0.3
        mesh_partition : int, optional
            mesh partition, by default 30
        strain : list, optional
            strain, by default [0.1, 0.0, 0.0]
        strain_amplitude : list, optional
            strain amplitude, by default None
        num_steps : int, optional
            number of steps, by default 100
        simulation_time : float, optional
            simulation time, by default 1.0
        num_cpu : int, optional
            number of cpu been used, by default 1
        hardening_law : Any, optional
            hardening law, by default LinearHardeningLaw()
        print_info : bool, optional
            print simulation information to the screen, by default False
        """

        # get the simulation information
        self.size = size
        self.radius = radius
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.strain_amplitude = strain_amplitude
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.hardening_law = hardening_law
        self.hardening_table = self.hardening_law.calculate_hardening_table()

        # update simulation information to logger
        self.logger.info("============== Simulation information ============")
        self.logger.info("size: {}".format(self.size))
        self.logger.info("radius: {}".format(radius))
        self.logger.info("young_modulus: {}".format(youngs_modulus))
        self.logger.info("poisson_ratio: {}".format(poisson_ratio))
        self.logger.info("mesh_partition: {}".format(mesh_partition))
        self.logger.info("strain: {}".format(strain))
        self.logger.info("strain_amplitude: {}".format(strain_amplitude))
        self.logger.info("hardening_law: {}".format(hardening_law))
        self.logger.info("hardening_table: {}".format(self.hardening_table))
        self.logger.info("num_steps: {}".format(num_steps))
        self.logger.info("simulation_time: {}".format(simulation_time))
        self.logger.info("num_cpu: {}".format(num_cpu))

        # have the simulation paras
        self.sim_paras = {"job_name": "hollow_plate",
                          "radius": radius,
                          "size": size,
                          "youngs_modulus": youngs_modulus,
                          "poisson_ratio": poisson_ratio,
                          "hardening_table": self.hardening_table,
                          "mesh_partition": mesh_partition,
                          "strain": strain,
                          "strain_amplitude": strain_amplitude,
                          "num_steps": num_steps,
                          "simulation_time": simulation_time,
                          "num_cpu": num_cpu}

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)
