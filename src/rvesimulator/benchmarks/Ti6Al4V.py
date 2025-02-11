#                                                                       Modules
# =============================================================================
# Standard
from .py3rve_base import Py3RVEBase
from rvesimulator.microstructure.sphere_particles import (
    SphereParticles, MicrostructureGenerator)
from rvesimulator.microstructure.circle_particles import CircleParticles
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict
import numpy as np
# local
import rvesimulator
from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
from rvesimulator.additions.hardening_law import (
    SwiftHardeningLaw, HardeningLaw)


#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


class Ti6Al4V_2D(Py3RVEBase):
    """Interface between python and abaqus of the CDDM RVE case

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
            "sim_script": "benchmark_abaqus_scripts.two_materials_rve",
            "sim_func": "VonMisesPlasticMatrixFiberRegularLoads",
            "post_script": "basic_analysis_scripts.post_process",
            "post_func": "PostProcess2D",
        }

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        youngs_modulus_matrix: float = 100.0,
        poisson_ratio_matrix: float = 0.3,
        youngs_modulus_fiber: float = 1.0,
        poisson_ratio_fiber: float = 0.19,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        hardening_law_fiber: HardeningLaw = SwiftHardeningLaw(),
        hardening_law_matrix: HardeningLaw = SwiftHardeningLaw(),
        seed: Any = None,
        print_info: bool = False,
        min_dist_factor: float = 1.2,
    ) -> None:
        """regular rve for Ti6Al4V

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
        hardening_law_fiber : Any, optional
            hardening law for the fiber material,
            by default LinearHardeningLaw()
        haderning_law_matrix : Any, optional
            hardening law for the matrix material,
            by default LinearHardeningLaw()
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get simulation information
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        self.youngs_modulus_matrix = youngs_modulus_matrix
        self.poisson_ratio_matrix = poisson_ratio_matrix
        self.youngs_modulus_fiber = youngs_modulus_fiber
        self.poisson_ratio_fiber = poisson_ratio_fiber
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.hardening_law_fiber = hardening_law_fiber
        self.hardening_law_matrix = hardening_law_matrix
        self.seed = seed
        # get hardening law
        self.hardening_table_fiber = \
            hardening_law_fiber.calculate_hardening_table()
        self.hardening_table_matrix = \
            hardening_law_matrix.calculate_hardening_table()

        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "youngs_modulus_matrix": youngs_modulus_matrix,
            "poisson_ratio_matrix": poisson_ratio_matrix,
            "youngs_modulus_fiber": youngs_modulus_fiber,
            "poisson_ratio_fiber": poisson_ratio_fiber,
            "hardening_table_fiber": self.hardening_table_fiber,
            "hardening_table_matrix": self.hardening_table_matrix,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu}

        self.mini_dist_factor = min_dist_factor

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        self.sim_info = {
            "job_name": "cddm_rve",
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
            "youngs_modulus_matrix": self.youngs_modulus_matrix,
            "poisson_ratio_matrix": self.poisson_ratio_matrix,
            "youngs_modulus_fiber": self.youngs_modulus_fiber,
            "poisson_ratio_fiber": self.poisson_ratio_fiber,
            "mesh_partition": self.mesh_partition,
            "hardening_table_fiber": self.hardening_table_fiber,
            "hardening_table_matrix": self.hardening_table_matrix,
            "num_steps": self.num_steps,
            "simulation_time": self.simulation_time,
            "strain": self.strain,
            "num_cpu": self.num_cpu,
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
        # create microstructure
        self.microstructure = CircleParticles(
            length=self.size,
            width=self.size,
            radius_mu=self.radius_mu,
            radius_std=self.radius_std,
            vol_req=self.vol_req,
            dist_min_factor=self.mini_dist_factor,
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


class Ti6Al4V_3D(Py3RVEBase):
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
            "sim_script": "benchmark_abaqus_scripts.j2_plastic_3d_rve",
            "sim_func": "j2_plastic_3d_rve",
            "post_script": "benchmark_abaqus_scripts.j2_plastic_3d_rve",
            "post_func": "post_process",
        }

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        youngs_modulus_matrix: float = 100.0,
        poisson_ratio_matrix: float = 0.3,
        youngs_modulus_fiber: float = 1.0,
        poisson_ratio_fiber: float = 0.19,
        mesh_partition: int = 100,
        strain: list = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 4,
        hardening_law_fiber: HardeningLaw = SwiftHardeningLaw(),
        hardening_law_matrix: HardeningLaw = SwiftHardeningLaw(),
        seed: Any = None,
        print_info: bool = False,
        min_dist_factor: float = 1.025,
        regular_load: bool = False,
        strain_amplitude: list = None,
        specify_microstructure: bool = False,
        microstructure_info: Dict = None,
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
        hardening_law_fiber : Any, optional
            hardening law for the fiber material,
            by default LinearHardeningLaw()
        haderning_law_matrix : Any, optional
            hardening law for the matrix material,
            by default LinearHardeningLaw()
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get simulation information
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        self.youngs_modulus_matrix = youngs_modulus_matrix
        self.poisson_ratio_matrix = poisson_ratio_matrix
        self.youngs_modulus_fiber = youngs_modulus_fiber
        self.poisson_ratio_fiber = poisson_ratio_fiber
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.hardening_law_fiber = hardening_law_fiber
        self.hardening_law_matrix = hardening_law_matrix
        self.seed = seed
        # specify microstructure
        self.specify_microstructure = specify_microstructure
        self.microstructure_info = microstructure_info
        # loading
        self.regular_load = regular_load
        self.strain_amplitude = strain_amplitude
        # get hardening law
        self.hardening_table_fiber = \
            hardening_law_fiber.calculate_hardening_table()
        self.hardening_table_matrix = \
            hardening_law_matrix.calculate_hardening_table()
        if self.regular_load:
            self.sim_paras = {
                "size": size,
                "radius_mu": radius_mu,
                "radius_std": radius_std,
                "vol_req": vol_req,
                "youngs_modulus_matrix": youngs_modulus_matrix,
                "poisson_ratio_matrix": poisson_ratio_matrix,
                "youngs_modulus_fiber": youngs_modulus_fiber,
                "poisson_ratio_fiber": poisson_ratio_fiber,
                "hardening_table_fiber": self.hardening_table_fiber,
                "hardening_table_matrix": self.hardening_table_matrix,
                "mesh_partition": mesh_partition,
                "strain": strain,
                "num_steps": num_steps,
                "simulation_time": simulation_time,
                "num_cpu": num_cpu}
        else:
            self.sim_paras = {
                "size": size,
                "radius_mu": radius_mu,
                "radius_std": radius_std,
                "vol_req": vol_req,
                "youngs_modulus_matrix": youngs_modulus_matrix,
                "poisson_ratio_matrix": poisson_ratio_matrix,
                "youngs_modulus_fiber": youngs_modulus_fiber,
                "poisson_ratio_fiber": poisson_ratio_fiber,
                "hardening_table_fiber": self.hardening_table_fiber,
                "hardening_table_matrix": self.hardening_table_matrix,
                "mesh_partition": mesh_partition,
                "strain": strain,
                "num_steps": num_steps,
                "simulation_time": simulation_time,
                "num_cpu": num_cpu,
                "strain_amplitude": strain_amplitude}

        self.mini_dist_factor = min_dist_factor

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        if self.regular_load:
            self.sim_info = {
                "job_name": "Ti6Al4V_3D",
                "location_information": self.microstructure_info[
                    "location_information"
                ],
                "radius_mu": self.microstructure_info["radius_mu"],
                "radius_std": self.microstructure_info[
                    "radius_std"],
                "len_start": self.microstructure_info["len_start"],
                "len_end": self.microstructure_info["len_end"],
                "wid_start": self.microstructure_info["wid_start"],
                "wid_end": self.microstructure_info["wid_end"],
                "hei_start": self.microstructure_info["hei_start"],
                "hei_end": self.microstructure_info["hei_end"],
                "youngs_modulus_matrix": self.youngs_modulus_matrix,
                "poisson_ratio_matrix": self.poisson_ratio_matrix,
                "youngs_modulus_fiber": self.youngs_modulus_fiber,
                "poisson_ratio_fiber": self.poisson_ratio_fiber,
                "mesh_partition": self.mesh_partition,
                "hardening_table_fiber": self.hardening_table_fiber,
                "hardening_table_matrix": self.hardening_table_matrix,
                "num_steps": self.num_steps,
                "simulation_time": self.simulation_time,
                "strain": self.strain,
                "num_cpu": self.num_cpu,
            }
        else:
            self.sim_info = {
                "job_name": "Ti6Al4V_3D",
                "location_information": self.microstructure_info[
                    "location_information"
                ],
                "radius_mu": self.microstructure_info["radius_mu"],
                "radius_std": self.microstructure_info[
                    "radius_std"],
                "len_start": self.microstructure_info["len_start"],
                "len_end": self.microstructure_info["len_end"],
                "wid_start": self.microstructure_info["wid_start"],
                "wid_end": self.microstructure_info["wid_end"],
                "hei_start": self.microstructure_info["hei_start"],
                "hei_end": self.microstructure_info["hei_end"],
                "youngs_modulus_matrix": self.youngs_modulus_matrix,
                "poisson_ratio_matrix": self.poisson_ratio_matrix,
                "youngs_modulus_fiber": self.youngs_modulus_fiber,
                "poisson_ratio_fiber": self.poisson_ratio_fiber,
                "mesh_partition": self.mesh_partition,
                "hardening_table_fiber": self.hardening_table_fiber,
                "hardening_table_matrix": self.hardening_table_matrix,
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

        if self.specify_microstructure:
            # plot microstructure
            MicrostructureGenerator().sphere_plot(
                fibers=np.array(
                    self.microstructure_info["location_information"]),
                length=self.microstructure_info["len_end"] -
                self.microstructure_info["len_start"],
                width=self.microstructure_info["wid_end"] -
                self.microstructure_info["wid_start"],
                height=self.microstructure_info["hei_end"] -
                self.microstructure_info["hei_start"],
                vol_frac=self.microstructure_info["vol_frac"],
                save_figure=True,
                fig_name="{}.png".format(self.microstructure_info["vol_frac"]))
            self.vol_frac = self.microstructure_info["vol_frac"]

        else:
            self.microstructure = SphereParticles(
                length=self.size,
                width=self.size,
                height=self.size,
                radius_mu=self.radius_mu,
                radius_std=self.radius_std,
                vol_req=self.vol_req,
                dist_min_factor=self.mini_dist_factor,
                num_cycle_max=20,
                num_guess_max=100000,
                num_fiber_max=20000,
            )
            # create microstructure
            self.microstructure.generate_microstructure(seed=self.seed)
            self.microstructure.to_abaqus_format()
            self.microstructure.plot_microstructure(save_figure=True,
                                                    fig_name="3d_rve_{}.png".
                                                    format(self.seed))
            self.vol_frac = self.microstructure.vol_frac
            self.microstructure_info = self.microstructure.microstructure_info
        self.logger.info("volume fraction: {}".format(self.vol_frac))
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


class Ti6Al4VElastic2D(Py3RVEBase):
    """Interface between python and abaqus of the 2D RVE case:
    Ti6Al4V plasticity for matrix material and elastic for fiber material

    Parameters
    ----------
    SimulationBase : class
        base class for simulation
    """

    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""

        logging.basicConfig(level=logging.INFO, filename="Ti6Al4VElastic.log")
        self.logger = logging.getLogger("abaqus_simulation")

        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "benchmark_abaqus_scripts.two_materials_rve",
            "sim_func": "VonMisesPlasticElasticPathLoads",
            "post_script": "basic_analysis_scripts.post_process",
            "post_func": "PostProcess2D",
        }

    def update_sim_info(
        self,
        size: float = 1.0,
        radius_mu: float = 0.10,
        radius_std: float = 0.01,
        vol_req: float = 0.30,
        youngs_modulus_matrix: float = 110000.0,
        poisson_ratio_matrix: float = 0.33,
        youngs_modulus_fiber: float = 22000.0,
        poisson_ratio_fiber: float = 0.19,
        mesh_partition: int = 200,
        strain: list = [0.05, 0.05, 0.05],
        strain_amplitude: list = None,
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 4,
        hardening_law: HardeningLaw = SwiftHardeningLaw(
            a=700, b=0.5, yield_stress=900),
        seed: Any = None,
        mini_dist_factor: float = 1.3,
        print_info: bool = False,
    ) -> None:
        """path dependent rve for cddm

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
        strain_amplitude : list, optional
            applied strain amplitude to mimic path dependence, by default None
        num_steps : int, optional
            number of simulation steps, by default 100
        simulation_time : float, optional
            total simulation time, by default 1.0
        num_cpu : int, optional
            number of cpu used for simulation, by default 1
        platform : str, optional
            platform for simulation, by default "ubuntu"
        hardening_law : class, optional
            hardening law for the simulation, by default LinearHardeningLaw()
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get simulation information
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        self.youngs_modulus_matrix = youngs_modulus_matrix
        self.poisson_ratio_matrix = poisson_ratio_matrix
        self.youngs_modulus_fiber = youngs_modulus_fiber
        self.poisson_ratio_fiber = poisson_ratio_fiber
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.strain_amplitude = strain_amplitude
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.hardening_law = hardening_law
        self.seed = seed
        self.mini_dist_factor = mini_dist_factor
        # get hardening law
        self.hardening_table = hardening_law.calculate_hardening_table()

        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "youngs_modulus_matrix": youngs_modulus_matrix,
            "poisson_ratio_matrix": poisson_ratio_matrix,
            "youngs_modulus_fiber": youngs_modulus_fiber,
            "poisson_ratio_fiber": poisson_ratio_fiber,
            "hardening_table": self.hardening_table,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "mini_dist_factor": mini_dist_factor, }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        self.sim_info = {
            "job_name": "cddm_rve",
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
            "youngs_modulus_matrix": self.youngs_modulus_matrix,
            "poisson_ratio_matrix": self.poisson_ratio_matrix,
            "youngs_modulus_fiber": self.youngs_modulus_fiber,
            "poisson_ratio_fiber": self.poisson_ratio_fiber,
            "mesh_partition": self.mesh_partition,
            "hardening_table": self.hardening_table,
            "num_steps": self.num_steps,
            "simulation_time": self.simulation_time,
            "strain": self.strain,
            "strain_amplitude": self.strain_amplitude,
            "num_cpu": self.num_cpu,
            "mini_dist_factor": self.mini_dist_factor,
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
        # create microstructure
        self.microstructure = CircleParticles(
            length=self.size,
            width=self.size,
            radius_mu=self.radius_mu,
            radius_std=self.radius_std,
            vol_req=self.vol_req,
            dist_min_factor=self.mini_dist_factor,
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
        # update logger on samples
        self.logger.info("==============        update info      ============")
        self.logger.info("sample: {}".format(sample))

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
