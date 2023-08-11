#                                                                       Modules
# =============================================================================
# Standard
import os

import rvesimulator
from rvesimulator.additions.hardening_law import HardeningLaw
from rvesimulator.additions.microstructure_wrapper import \
    CircleSVEMicroStructure

from .shared_functionalities import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================

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
            "sim_path": "benchmark_abaqus_scripts.two_materials_sve",
            "sim_script": "VonMisesPlasticElasticRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.15,
        sve_geometry_benchmark: str = "task1",
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
        print_info: bool = False,
        **kwargs,
    ) -> None:
        # get the micro_structure information
        micro_structure_generator = CircleSVEMicroStructure(
            size=size,
            radius=radius,
        )
        microstructure_info = micro_structure_generator.location_information(
            task=sve_geometry_benchmark
        )
        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardening_table = law_function(**kwargs)

        self.sim_paras = {"job_name": "two_materials_sve",
                          "size": size,
                          "location_information": microstructure_info,
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
            "job_name": "two_materials_sve",
            "size": size,
            "location_information": microstructure_info,
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
            "sim_path": "benchmark_abaqus_scripts.two_materials_sve",
            "sim_script": "VonMisesPlasticElasticPathLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }

        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.15,
        sve_geometry_benchmark: str = "task1",
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
        print_info: bool = False,
        **kwargs,
    ) -> None:
        # get the micro_structure information
        micro_structure_generator = CircleSVEMicroStructure(
            size=size,
            radius=radius,
        )
        microstructure_info = micro_structure_generator.location_information(
            task=sve_geometry_benchmark
        )
        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardening_table = law_function(**kwargs)

        self.sim_paras = {"job_name": "two_materials_sve",
                          "location_information": microstructure_info,
                          "size": size,
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
            "job_name": "two_materials_sve",
            "location_information": microstructure_info,
            "size": size,
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
