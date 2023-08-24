import json
import os

from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator

# local functions
from .utils import create_dir


class SimulationBase:
    """base class of rve simulation problems"""

    def run_simulation(
        self,
        sample: dict = None,
        folder_index: int = None,
        sub_folder_index: int = None,
        third_folder_index: int = None,
    ) -> dict:
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
            sub_folder_index,
            third_folder_index,
        )
        # update the geometry info for microstructure
        self._update_sample_info(sample=sample)
        # change folder to main folder
        os.chdir(self.main_folder)
        simulator = AbaqusSimulator(
            sim_info=self.sim_info, folder_info=self.folder_info
        )
        # run abaqus simulation
        simulator.run()
        # get the simulation results back
        results = simulator.read_back_results()

        return results

    def run_batch_simulation(self) -> any:
        """run batch simulation

        Returns
        -------
        any
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
