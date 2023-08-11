#                                                                       Modules
# =============================================================================
# Standard
import os
from typing import Any

import rvesimulator
from rvesimulator.additions.hardening_law import HardeningLaw
from rvesimulator.additions.microstructure_wrapper import CircleMircoStructure

from .shared_functionalities import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


class VonMisesPlasticElasticRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.two_materials_rve",
            "sim_script": "VonMisesPlasticElasticRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

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
        platform: str = "ubuntu",
        hardening_law: str = "swift",
        seed: any = None,
        print_info: bool = False,
        **kwargs,
    ) -> None:
        """update simulation information for CDDM case

        Parameters
        ----------
        size : float, optional
            size of rve, by default 0.048
        radius_mu : float, optional
            radius mean of rve, by default 0.003
        radius_std : float, optional
            radius standard deviation of rve, by default 0.0
        vol_req : float, optional
            volume requirement of rve, by default 0.30
        youngs_modulus_matrix : float, optional
            youngs modulus of matrix material, by default 100.0
        poisson_ratio_matrix : float, optional
            poisson ratio of matrix material, by default 0.3
        youngs_modulus_fiber : float, optional
            youngs modulus of fiber material, by default 1.0
        poisson_ratio_fiber : float, optional
            poisson ratio of fiber material, by default 0.19
        mesh_partition : int, optional
            mesh partitions for edges, by default 30
        strain : list, optional
            applied maximum strain, by default [0.1, 0.0, 0.0]
        num_steps : int, optional
            number of simulation steps, by default 100
        simulation_time : float, optional
            total simulation time, by default 1.0
        num_cpu : int, optional
            number of cpu been used, by default 1
        platform : str, optional
            simulation platform, by default "ubuntu"
        hardening_law : str, optional
            name of hardening law, by default "swift"
        seed : any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )
        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardening_table = law_function(**kwargs)

        self.sim_paras = {"job_name": "two_materials_rve",
                          "size": size,
                          "radius_mu": radius_mu,
                          "radius_std": radius_std,
                          "vol_req": vol_req,
                          "vol_frac": vol_frac,
                          "youngs_modulus_matrix": youngs_modulus_matrix,
                          "poisson_ratio_matrix": poisson_ratio_matrix,
                          "hardening_law": hardening_law,
                          "hardening_paras": kwargs,
                          "youngs_modulus_fiber": youngs_modulus_fiber,
                          "poisson_ratio_fiber": poisson_ratio_fiber,
                          "mesh_partition": mesh_partition,
                          "strain": strain,
                          "num_steps": num_steps,
                          "simulation_time": simulation_time,
                          "num_cpu": num_cpu,
                          "platform": platform, }

        self.sim_info = {
            "job_name": "two_materials_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "youngs_modulus_matrix": youngs_modulus_matrix,
            "poisson_ratio_matrix": poisson_ratio_matrix,
            "youngs_modulus_fiber": youngs_modulus_fiber,
            "poisson_ratio_fiber": poisson_ratio_fiber,
            "mesh_partition": mesh_partition,
            "hardening_table": hardening_table,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

# =============================================================================


class VonMisesPlasticElasticPathLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.two_materials_rve",
            "sim_script": "VonMisesPlasticElasticPathLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

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
        strain_amplitude: list = None,
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        hardening_law: str = "swift",
        seed: Any = None,
        print_info: bool = False,
        **kwargs,
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
        hardening_law : str, optional
            name of hardening law, by default "swift"
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )
        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardening_table = law_function(**kwargs)

        # simulation paras
        self.sim_paras = {"job_name": "two_materials_rve",
                          "size": size,
                          "radius_mu": radius_mu,
                          "radius_std": radius_std,
                          "vol_req": vol_req,
                          "vol_frac": vol_frac,
                          "youngs_modulus_matrix": youngs_modulus_matrix,
                          "poisson_ratio_matrix": poisson_ratio_matrix,
                          "hardening_law": hardening_law,
                          "hardening_paras": kwargs,
                          "youngs_modulus_fiber": youngs_modulus_fiber,
                          "poisson_ratio_fiber": poisson_ratio_fiber,
                          "mesh_partition": mesh_partition,
                          "strain": strain,
                          "strain_amplitude": strain_amplitude,
                          "num_steps": num_steps,
                          "simulation_time": simulation_time,
                          "num_cpu": num_cpu,
                          "platform": platform, }

        self.sim_info = {
            "job_name": "two_materials_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "youngs_modulus_matrix": youngs_modulus_matrix,
            "poisson_ratio_matrix": poisson_ratio_matrix,
            "youngs_modulus_fiber": youngs_modulus_fiber,
            "poisson_ratio_fiber": poisson_ratio_fiber,
            "mesh_partition": mesh_partition,
            "hardening_table": hardening_table,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "strain_amplitude": strain_amplitude,
            "num_cpu": num_cpu,
            "platform": platform,
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)
