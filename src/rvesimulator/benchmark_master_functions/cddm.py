#                                                                       Modules
# =============================================================================
# Standard
import json
import os

import rvesimulator
from rvesimulator.additions.hardening_law import HardeningLaw
from rvesimulator.additions.microstructure_wrapper import (
    CircleMircoStructure,
)
from .shared_functionalities import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


# =============================================================================
# class HardeningLawInterface(Protocol):
#     def __init__(self, *args, **kwargs):
#         ...


#     def to_dict(self):  
#         return {'table': self.table} # self.__dict__


# class LinearHardeningLaw(HardeningLawInterface):
#     def __init__(self,a, b):

#         ...

#         self.table 
                     
    


# class SimulatorInfo:
#     def __init__(self, folder_info: FolderInfo, hardening_law: HardeningLawInterface = None, microstructure: MicrostructureInterface):
#         self.folder_info = folder_info
#         self.hardening_law = hardening_law
#         self.microstructure = microstructure

#         self.run_checks()

#     def run_checks(self):
#         ...

#     def to_dict(self) -> dict:
#         sim_info: dict = {}

#         sim_info.update(self.hardening_law.to_dict())
#         sim_info.update(self.microstructure.to_dict())

#         ...

#         return sim_info

# class VonMisesPlasticElasticRegularLoads2(SimulatorInfo):
#     def run_checks(self):
#         # Check if the arguments are valid
#         assert isinstance(self.hardening_law, LinearHardeningLaw)

# ###

# sim_info = SimulatorInfo.to_dict()
# folder_info = FolderInfo.to_dict()
# class AbaqusSimulator:
#     def __init__(self, simulator_info = SimulatorInfo, FolderInfo):
#         ...

#     def execute



class VonMisesPlasticElasticRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + "/scriptbase",
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
        poisson_ratio_fiber: float  = 0.19, 
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
            "script_path": os.path.dirname(rvesimulator.__file__) + "/scriptbase",
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
        poisson_ratio_fiber: float  = 0.19, 
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        strain_amplitude: list = None,
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        hardening_law: str = "swift",
        seed: any = None,
        print_info: bool = False,
        **kwargs,
    ) -> None:
        
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
                          "strain_amplitude":strain_amplitude,
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
