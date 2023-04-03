import os

import numpy as np
from f3dasm.design import ExperimentData
from f3dasm.simulation.abaqus_simulator import AbaqusSimulator

# local functions
from .utils import create_dir


class SimulationBase:
    """base class of rve simulation problems"""

    def __init__(self) -> None:
        """initialization"""
        pass

    def run_f3dasm(self, data: ExperimentData) -> ExperimentData:
        """run simulation via f3dasm manners

        Parameters
        ----------
        data : ExperimentData
            data object

        Returns
        -------
        ExperimentData
            updated data object
        """
        # get the samples
        samples = data.data.input.to_dict("record")
        for ii in range(len(data.data)):
            results = self.run_simulation(
                sample=samples[ii], third_folder_index=ii
            )
            # fill the data class
            data.data["output"] = data.data["output"].astype(object)
            for jj in range(len(list(data.data["output"].keys()))):
                data.data[("output", list(data.data["output"].keys())[jj])][
                    ii
                ] = results[list(data.data["output"].keys())[jj]]

        return data

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
            a dict contains the information of design variables, by default None
        folder_index : int, optional
            first fodler index, by default None
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
        """create folders for excuting abaqus simulations

        Parameters
        ----------
        folder_index : int, optional
            first folder index , by default None
        sub_folder_index : _type_, optional
            second folder index , by default None
        third_folder_index : _type_, optional
            third forlder index, by default None

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
                raise ValueError("provide sub_folder_index")
            elif third_folder_index is None:
                raise ValueError("provide third_folder_index")
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
