import json
import os
from typing import Any

# local functions
from .utils import create_dir


class SimulationBase:
    """base class of rve simulation problems"""

    def run_simulation(self) -> Any:
        """run simulation

        Returns
        -------
        Any
            should be defined according to different problems
        """
        raise NotImplementedError("should be implemented in subclass")

    def run_batch_simulation(self) -> Any:
        """run batch simulation

        Returns
        -------
        Any
            should be defined according to different problems

        Raises
        ------
        NotImplementedError
            should be implemented in specific problems
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
        folder_index: int = None,
        sub_folder_index: int = None,
        third_folder_index: int = None,
    ) -> None:
        """create folders for executing abaqus simulations, it's a third
        layer folder structure, you have to provide at least one folder index

        Parameters
        ----------
        folder_index : int, optional
            first folder index , by default None
        sub_folder_index : int, optional
            second folder index , by default None
        third_folder_index : int, optional
            third folder index, by default None

        Raises
        ------
        ValueError
            provide sub_folder_index
        ValueError
            provide third_folder_index
        """
        if folder_index is None:
            if sub_folder_index is None:
                self.folder_info["current_work_directory"] = "case_" + str(
                    third_folder_index
                )
            else:
                if third_folder_index is None:
                    self.folder_info[
                        "current_work_directory"
                    ] = "point_" + str(sub_folder_index)
                else:
                    self.folder_info["current_work_directory"] = (
                        "point_"
                        + str(sub_folder_index)
                        + "/case_"
                        + str(third_folder_index)
                    )
        else:
            if sub_folder_index is None:
                if third_folder_index is None:
                    self.folder_info[
                        "current_work_directory"
                    ] = "point_" + str(folder_index)
                else:
                    self.folder_info["current_work_directory"] = (
                        "gen_" + str(folder_index)
                        + "/case_" + str(third_folder_index))
            elif third_folder_index is None:
                self.folder_info["current_work_directory"] = (
                    "gen_" + str(folder_index)
                    + "/point_" + str(sub_folder_index))
            else:
                self.folder_info["current_work_directory"] = (
                    "gen_"
                    + str(folder_index)
                    + "/point_"
                    + str(sub_folder_index)
                    + "/case_"
                    + str(third_folder_index)
                )
        new_path = create_dir(
            current_folder=self.folder_info["main_work_directory"],
            dir_name=self.folder_info["current_work_directory"],
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
