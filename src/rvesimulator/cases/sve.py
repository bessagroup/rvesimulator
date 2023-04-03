#                                                                       Modules
# =============================================================================
# Standard
import json
import os

import rvesimulator
from rvesimulator.additionals.hardening_law import HardeningLaw
from rvesimulator.additionals.microstructure_wrapper import (
    CircleSVEMicroStructure,
)

from .base import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class SVE(SimulationBase):
    def __init__(self) -> None:
        self.main_folder = os.getcwd()
        self.samples = None
        self.rve_geometry = None
        self.abaqus_paras = None
        self.sim_info = None
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__),
            "current_work_directory": "point_1",
            "sim_path": "scriptbase.composite_sve_script",
            "sim_script": "SVE",
            "post_path": "scriptbase.basic_script.postprocess",
            "post_script": "RVEPostProcess2D",
        }
        self.vol_req = None
        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.10,
        task: str = "task1",
        mesh_partition: int = 100,
        strain: list = [0.05, 0.02, 0.02],
        strain_amplitude: list = None,
        simulation_time: float = 1.0,
        num_steps: int = 100,
        E_matrix: float = 100.0,
        Pr_matrix: float = 0.3,
        E_fiber: float = 1.0,
        Pr_fiber: float = 0.19,
        hardening_law: str = "swift",
        num_cpu: int = 1,
        platform: str = "ubuntu",
        print_info: bool = False,
        **kwargs,
    ) -> None:
        """update simulation information

        Parameters
        ----------
        size : float, optional
            size of rve, by default 0.048
        radius : float, optional
            radiusof fibers, by default 0.003
        task: str
            flag to get microstructure information
        mesh_partition : int, optional
            mesh partition, by default 100
        strain : list, optional
            applied strain, by default [0.05, 0.0, 0.0]
        strain_amplitude : list, optional
            strain amplitude, by default None
        simulation_time : float, optional
            simulation time, by default 1.0
        num_steps : int, optional
            number steps , by default 100
        E_matrix : float, optional
            young's modulus of matrix material , by default 100.0
        Pr_matrix : float, optional
            poisson ratio of matrix material, by default 0.3
        E_fiber : float, optional
            young's modulsu of fiber material, by default 1.0
        Pr_fiber : float, optional
            poisson ratio of fiber material, by default 0.19
        hardening_law : str, optional
            name of hardening law, by default "swift"
        num_cpu : int, optional
            number of cpu, by default 1
        platform : str, optional
            excution platform, by default "ubuntu"
        seed : any, optional
            seed of microstructure generation, by default None
        print_info : bool, optional
            print info, by default False
        """
        # generate the microstructure
        mircostructure_generator = CircleSVEMicroStructure(
            size=size,
            radius=radius,
        )

        microstructure_info = mircostructure_generator.location_information(
            task=task
        )

        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardenning_table = law_function(**kwargs)

        self.sim_paras = {
            "length": size,
            "width": size,
            "radius": radius,
            "task": task,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "strain_amplitude": strain_amplitude,
            "simulation_time": simulation_time,
            "num_steps": num_steps,
            "num_cpu": num_cpu,
            "platform": platform,
            "hardening_law": hardening_law,
            "hardening paras": kwargs,
            "E_matrix": E_matrix,
            "Pr_matrix": Pr_matrix,
            "E_fiber": E_fiber,
            "Pr_fiber": Pr_fiber,
        }

        if print_info is True:
            print("Simulation information: \n")
            print(json.dumps(self.sim_paras, indent=4))

        self.sim_info = {
            "length": size,
            "width": size,
            "location_information": microstructure_info,
            "E_matrix": E_matrix,
            "Pr_matrix": Pr_matrix,
            "E_fiber": E_fiber,
            "Pr_fiber": Pr_fiber,
            "hardening_table": hardenning_table,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "strain_amplitude": strain_amplitude,
            "simulation_time": simulation_time,
            "num_steps": num_steps,
            "job_name": "sve",
            "num_cpu": num_cpu,
            "platform": platform,
        }
