"""here"""
#                                                                       Modules
# =============================================================================
# Standard
import json
import os

import rvesimulator
from rvesimulator.additions.hardening_law import HardeningLaw

from .shared_functionalities import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


# =============================================================================


class ElasticRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""

        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.hollow_plate_sve",
            "sim_script": "ElasticRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(print_info=False)

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        num_cpu: int = 1,
        platform: str = "ubuntu",
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
        self.sim_info = {
            "job_name": "hollow_plate",
            "radius": radius,
            "size": size,
            "youngs_modulus": youngs_modulus,
            "poisson_ratio": poisson_ratio,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_info)

# =============================================================================
class VonMisesPlasticRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.hollow_plate_sve",
            "sim_script": "VonMisesPlasticRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        hardening_law: str = "swift",
        print_info: bool = False,
        **kwargs,
    ) -> None:
        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardening_table = law_function(**kwargs)

        self.sim_paras = {"job_name": "hollow_plate",
                          "radius": radius,
                          "size": size,
                          "youngs_modulus": youngs_modulus,
                          "poisson_ratio": poisson_ratio,
                          "hardening_law": hardening_law,
                          "hardening_paras": kwargs,
                          "mesh_partition": mesh_partition,
                          "strain": strain,
                          "num_steps": num_steps,
                          "simulation_time": simulation_time,
                          "num_cpu": num_cpu,
                          "platform": platform, }

        self.sim_info = {
            "job_name": "hollow_plate",
            "radius": radius,
            "size": size,
            "youngs_modulus": youngs_modulus,
            "poisson_ratio": poisson_ratio,
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
class VonMisesPlasticPathLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.hollow_plate_sve",
            "sim_script": "VonMisesPlasticPathLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

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
        platform: str = "ubuntu",
        hardening_law: str = "swift",
        print_info: bool = False,
        **kwargs,
    ) -> None:
        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardening_table = law_function(**kwargs)

        self.sim_paras = {"job_name": "hollow_plate",
                          "radius": radius,
                          "size": size,
                          "youngs_modulus": youngs_modulus,
                          "poisson_ratio": poisson_ratio,
                          "hardening_law": hardening_law,
                          "hardening_paras": kwargs,
                          "mesh_partition": mesh_partition,
                          "strain": strain,
                          "strain_amplitude": strain_amplitude,
                          "num_steps": num_steps,
                          "simulation_time": simulation_time,
                          "num_cpu": num_cpu,
                          "platform": platform, }

        self.sim_info = {
            "job_name": "hollow_plate",
            "radius": radius,
            "size": size,
            "youngs_modulus": youngs_modulus,
            "poisson_ratio": poisson_ratio,
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
