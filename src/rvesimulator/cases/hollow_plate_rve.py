#                                                                       Modules
# =============================================================================
# Standard
import json
import os

import rvesimulator

from .base import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class NaiveHollowPlate(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollowplate case"""

        self.main_folder = os.getcwd()
        self.samples = None
        self.sim_info = None
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__),
            "current_work_directory": "point_1",
            "sim_path": "scriptbase.hollow_plate_script",
            "sim_script": "NaiveHollowPlate",
            "post_path": "scriptbase.basic_script.postprocess",
            "post_script": "RVEPostProcess2D",
        }

        self.update_sim_info(print_info=False)

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poisson_ratio: float = 0.3,
        mesh_portion: int = 30,
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
            platform of excution, by default "ubuntu"
        print_info : bool, optional
            print the simulation information, by default False
        """
        self.sim_info = {
            "job_name": "hollowplate",
            "radius": radius,
            "size": size,
            "youngs_modulus": youngs_modulus,
            "poisson_ratio": poisson_ratio,
            "mesh_portion": mesh_portion,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
        }
        if print_info is True:
            print("Simulation information: \n")
            print(json.dumps(self.sim_info, indent=4))
