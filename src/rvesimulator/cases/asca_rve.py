#                                                                       Modules
# =============================================================================
# Standard
import json
import os

import rvesimulator
from rvesimulator.additionals.hardening_law import HardeningLaw
from rvesimulator.additionals.microstructure_wrapper import (
    CircleMircoStructure,
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


class ASCARVE(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the asca rve case"""
        self.main_folder = os.getcwd()
        self.samples = None
        self.rve_geometry = None
        self.abaqus_paras = None
        self.sim_info = None
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__),
            "current_work_directory": None,
            "sim_path": "scriptbase.asca_rve_script",
            "sim_script": "ASCARVE",
            "post_path": "scriptbase.basic_script.postprocess",
            "post_script": "RVEPostProcess2D",
        }
        self.vol_req = None
        self.update_sim_info(yield_stress=0.5, a=0.2, b=0.4)

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.000,
        vol_req: float = 0.30,
        mesh_partition: int = 100,
        strain: list = [0.05, 0.0, 0.0],
        simulation_time: float = 1.0,
        num_steps: int = 100,
        E_matrix: float = 100.0,
        Pr_matrix: float = 0.3,
        E_fiber: float = 1.0,
        Pr_fiber: float = 0.19,
        hardening_law: str = "swift",
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: any = None,
        print_info: bool = False,
        **kwargs,
    ) -> None:
        """_summary_

        Parameters
        ----------
        size : float, optional
            size of rve, by default 0.048
        radius_mu : float, optional
            radius mean of fibers, by default 0.003
        radius_std : float, optional
            radius std of fibers , by default 0.000
        vol_req : float, optional
            required volume fraction, by default 0.30
        mesh_partition : int, optional
            mesh partition, by default 100
        strain : list, optional
            applied strain, by default [0.05, 0.0, 0.0]
        simulation_time : float, optional
            simulation time, by default 1.0
        num_steps : int, optional
            number steps , by default 100
        E_matrix : float, optional
            youngs modulus of matrix material, by default 100.0
        Pr_matrix : float, optional
            poisson ratio of matrix material, by default 0.3
        E_fiber : float, optional
            youngs modulus of fiber material, by default 1.0
        Pr_fiber : float, optional
            poisson ratio of fiber material, by default 0.19
        hardening_law : str, optional
            name of hardening law, by default "swift"
        num_cpu : int, optional
            cpu number, by default 1
        platform : str, optional
            platform of excution, by default "ubuntu"
        seed : any, optional
            seed for generating microstructure, by default None
        print_info : bool, optional
            print info, by default False
        """
        self.vol_req = vol_req
        # generate the microstructure
        mircostructure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mircostructure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )

        # generate the hardening law
        hardening = HardeningLaw()
        law_function = getattr(hardening, hardening_law)
        hardenning_table = law_function(**kwargs)

        self.sim_paras = {
            "length": size,
            "width": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "mesh_partition": mesh_partition,
            "strain": strain,
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
            "job_name": "asca_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "mesh_partition": mesh_partition,
            "E_matrix": E_matrix,
            "Pr_matrix": Pr_matrix,
            "E_fiber": E_fiber,
            "Pr_fiber": Pr_fiber,
            "hardening_table": hardenning_table,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
            "simulation_time": simulation_time,
            "num_steps": num_steps,
        }
