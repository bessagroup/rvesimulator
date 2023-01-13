#                                                                       Modules
# =============================================================================
# Standard
import os
import pickle

# local
import rvesimulator
from rvesimulator.simulators.abaqus_simulator import AbaqusSimulator
from rvesimulator.simulators.utils import create_dir

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class HollowPlateRVE:
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollowplate case"""

        self.main_folder = os.getcwd()
        self.samples = None
        self.sim_info = None
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__),
            "current_work_directory": "point_1",
            "sim_path": "scriptbase.hollow_plate_rve",
            "sim_script": "HollowPlateRVE",
            "post_path": "scriptbase.postprocess",
            "post_script": "RVEPostProcess",
        }

        self.update_sim_info(print_info=True)

    def run_simulation(self, data: dict) -> dict:
        """run simulation sequentially

        Parameters
        ----------
        data : dict
            design of experiments without responses

        Returns
        -------
        pandas.DataFrame
            design of experiments with responses
        """
        # number of samples
        self.data = data
        samples = data["samples"].to_dict("records")
        responses = data["responses"]
        num_samples = len(samples)
        for ii in range(num_samples):
            # update the simulation information
            self.folder_info["current_work_directory"] = "point_" + str(ii)
            new_path = create_dir(
                current_folder=self.folder_info["main_work_directory"],
                dirname=self.folder_info["current_work_directory"],
            )
            os.chdir(new_path)
            self._update_sim_info(samples[ii])
            os.chdir(self.main_folder)
            abaqus_wrapper = AbaqusSimulator(
                sim_info=self.sim_info, folder_info=self.folder_info
            )
            # run abaqus simulation
            abaqus_wrapper.run()
            # get the simulation results back
            results = abaqus_wrapper.read_back_results()
            # update DoE information
            for jj in range(len(list(responses.keys()))):
                responses.at[ii, list(responses.keys())[jj]] = results[
                    list(responses.keys())[jj]
                ]

            # save results
            self._save_data(responses)

        return self.data

    def update_sim_info(
        self,
        size: float = 1.0,
        radius: float = 0.2,
        youngs_modulus: float = 100.0,
        poission_ratio: float = 0.3,
        mesh_portion: int = 30,
        loads: list = [0.1, 0.0, 0.0],
        num_cpu: int = 1,
        platform: str = "ubuntu",
        print_info: bool = False,
    ) -> None:
        """update parameters

        Parameters
        ----------
        size : float, optional
            size of rve, by default 1.0
        radius : float, optional
            radius of the hole, by default 0.2
        youngs_modulus : float, optional
            Young's modulus , by default 100.0
        poission_ratio : float, optional
            Poission's ratio , by default 0.3
        mesh_portion : int, optional
            mesh portion, by default 30
        loads : list, optional
            loadings, by default [0.1, 0.0, 0.0]
        print_info : bool, optional
            print simulation inforation to screen, by default False
        """

        self.sim_info = {
            "job_name": "HollowPlateRVE",
            "radius": radius,
            "size": size,
            "youngs_modulus": youngs_modulus,
            "poission_ratio": poission_ratio,
            "mesh_portion": mesh_portion,
            "loads": loads,
            "num_cpu": num_cpu,
            "platform": platform,
        }
        if print_info is True:
            print(f"The simulation information is : {self.sim_info}")

    def _update_sim_info(self, sample) -> None:
        """update the design variables"""
        self.sim_info.update(sample)

    def _save_data(self, responses) -> None:
        """save data to json file"""
        self.data["responses"] = responses
        working_folder = os.getcwd()
        with open("data.pickle", "wb") as file:
            pickle.dump(self.data, file)
        os.chdir(working_folder)
