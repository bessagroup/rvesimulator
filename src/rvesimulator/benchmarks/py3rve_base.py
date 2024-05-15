"""
Module for simulation base class
"""
#                                                                       Modules
# =============================================================================
# Standard
import json
import os
from abc import ABC
from typing import Any, Dict

# local
from ..additions.utils import create_dir

#                                                          Authorship & Credits
# =============================================================================
__author__ = 'Jiaxiang Yi (J.Yi@tudelft.nl)'
__credits__ = ['Jiaxiang Yi']
__status__ = 'Stable'
# =============================================================================
#
# =============================================================================


class Py3RVEBase(ABC):
    """base class of rve simulation problems"""

    def __init__(self,
                 sim_info: Dict,
                 folder_info: Dict) -> None:

        self.sim_info = sim_info
        self.folder_info = folder_info

    def run_simulation(self) -> Any:
        """run simulation

        Returns
        -------
        Any
            should be defined according to different problems
        """
        raise NotImplementedError("should be implemented in subclass")

    def update_sim_info(self) -> None:
        """update the required simulation information

        Raises
        ------
        NotImplementedError
            should be implemented in specific problems
        """
        raise NotImplementedError("should be implemented in subclass")

    def _update_sample_info(self, sample: dict) -> None:
        """update the sample information

        Parameters
        ----------
        sample : dict
            sample
        """
        self.sim_info.update(sample)

    def _create_working_folder(
        self,
        folder_index: int,
    ) -> None:
        """create folders for executing abaqus simulations, it's a third
        layer folder structure, you have to provide at least one folder index

        """

        # update the current folder index
        self.folder_info["current_dir"] = "point_" + str(folder_index)

        # create a new folder
        new_path = create_dir(
            current_folder=self.folder_info["main_dir"],
            dir_name=self.folder_info["current_dir"],
        )
        self.working_folder = new_path

        os.chdir(new_path)

    def _print_sim_info(self, info: dict) -> None:
        """print simulation information to screen

        Parameters
        ----------
        info : dict
            a dict contains simulation information
        """

        print("Simulation information: \n")
        print(json.dumps(info, indent=4))
