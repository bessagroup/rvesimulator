#                                                                       Modules
# =============================================================================
# Standard
import json
import os
import pickle

# local
import rvesimulator
from rvesimulator.microstructures.circle_particles import (
    CircleParticles,
)
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


class AscaRVE:
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
            "current_work_directory": "point_1",
            "sim_path": "scriptbase.asca_rve",
            "sim_script": "ASCARVE",
            "post_path": "scriptbase.postprocess",
            "post_script": "RVEPostProcess2D",
        }
        self.vol_req = None
        self.update_sim_info()

    def run_simulation(
        self,
        data: dict = None,
        save_source_files: bool = True,
    ) -> dict:
        """run simulation sequentially

        Parameters
        ----------
        data : dict, optional
            design of experiments without responses
        save_source_files: bool, optional
            if True, it will create a new folder for every sample and
            more memory is needed for the simulation task
            other wise the simulation will be excuted in the Data/data
            folder and the simulation result will be save as Data/data.pickle
            file.

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
            if save_source_files is True:
                self.folder_info["current_work_directory"] = "point_" + str(ii)
                new_path = create_dir(
                    current_folder=self.folder_info["main_work_directory"],
                    dirname=self.folder_info["current_work_directory"],
                )
                os.chdir(new_path)
                log_file = "results.p"
                if os.path.exists(log_file):
                    print("remove results succesfully \n")
                    os.remove(log_file)
            else:
                self.folder_info["current_work_directory"] = "data"
                new_path = create_dir(
                    current_folder=self.folder_info["main_work_directory"],
                    dirname=self.folder_info["current_work_directory"],
                )
                os.chdir(new_path)
                log_file = "results.p"
                if os.path.exists(log_file):
                    print("remove results succesfully \n")
                    os.remove(log_file)
            # update the geometry info for microstructure
            self._update_sample_info(sample=samples[ii])

            # generating the micro-structure
            if os.path.isfile("micro_structure_info.json"):
                print("micro-structure file already exist\n")

            else:
                self.vol_frac = self.micro_structure_generation(
                    length=self.rve_geometry["length"],
                    width=self.rve_geometry["width"],
                    radius_mu=self.rve_geometry["radius_mu"],
                    radius_std=self.rve_geometry["radius_std"],
                    vol_req=self.vol_req,
                )
                print("micro-structure be generated successfully\n")

            list_iter = list(responses.keys())
            if "vol_frac" in list_iter:
                responses.at[ii, "vol_frac"] = self.vol_frac
                list_iter.remove("vol_frac")

            # complete simulation info for Abaqus
            self._complete_information()
            os.chdir(self.main_folder)

            # initilize the abaqus simulator
            abaqus_wrapper = AbaqusSimulator(
                sim_info=self.sim_info, folder_info=self.folder_info
            )
            abaqus_wrapper.run()
            results = abaqus_wrapper.read_back_results()
            # update DoE information
            for jj in range(len(list_iter)):
                responses.at[ii, list_iter[jj]] = results[list_iter[jj]]
            self._save_data(responses)

        return self.data

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0001,
        vol_req: float = 0.30,
        mesh_partition: int = 30,
        loads: list = [0.05, 0.0, 0.0],
        simulation_time: float = 1.0,
        print_info: bool = False,
        num_cpu: int = 1,
        platform: str = "ubuntu",
    ) -> None:
        """update parameters of asca rve simulation

        Parameters
        ----------
        size : float, optional
            the size of the RVE, by default 0.048
        radius : float, optional
            radius of , by default 0.003
        vol_req : float, optional
            volume requirement, by default 0.30
        mesh_partition : int, optional
            number of partition of every edge of the RVE, by default 30
        loads : list, optional
            a list of loads [Exx, Eyy, Exy], by default [0.05, 0.0, 0.0]
        simulation_time : float, optional
            total simulation time of Abaqus, by default 10.0
        print_info : bool, optional
            a flag to print info or not, by default False
        """

        self.vol_req = vol_req
        self.rve_geometry = {
            "length": size,
            "width": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
        }
        self.abaqus_paras = {
            "mesh_partition": mesh_partition,
            "loads": loads,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform,
        }
        if print_info:
            print(f"geometry information of RVE: {self.rve_geometry}")
            print(f"vol_req is: {self.vol_req}")
            print(f"Info of Abaqus simulation : {self.abaqus_paras} \n")

    def _update_sample_info(self, sample: dict) -> None:

        # update 'vol_req'
        if "vol_req" in sample.keys():
            self.vol_req = sample["vol_req"]

        # update 'size'
        if "size" in sample.keys():
            self.rve_geometry["length"] = sample["size"]
            self.rve_geometry["width"] = sample["size"]

        # update 'radius'
        if "radius" in sample.keys():
            self.rve_geometry["radius"] = sample["radius"]

        if "mesh_partition" in sample.keys():
            self.abaqus_paras["mesh_partition"] = sample["mesh_partition"]

        if "num_cpu" in sample.keys():
            self.abaqus_paras["num_cpu"] = sample["num_cpu"]

        if "platform" in sample.keys():
            self.abaqus_paras["platform"] = sample["platform"]

    def _complete_information(self) -> None:
        """
        This function is used to complete information for abaqus simulation

        Returns
        -------
        """

        # open the json file for micro-structures
        file = "micro_structure_info.json"
        with open(file, "r") as f:
            location_info = json.load(f)

        self.sim_info = {
            "job_name": "asca_rve",
            "location_information": location_info["location_information"],
            "radius_mu": location_info["radius_mu"],
            "radius_std": location_info["radius_std"],
            "len_start": location_info["len_start"],
            "len_end": location_info["len_end"],
            "wid_start": location_info["wid_start"],
            "wid_end": location_info["wid_end"],
            "mesh_partition": self.abaqus_paras["mesh_partition"],
            "loads": self.abaqus_paras["loads"],
            "simulation_time": self.abaqus_paras["simulation_time"],
            "num_cpu": self.abaqus_paras["num_cpu"],
            "platform": self.abaqus_paras["platform"],
        }

    @staticmethod
    def micro_structure_generation(
        length: float,
        width: float,
        radius_mu: float,
        radius_std: float,
        vol_req: float,
    ) -> float:
        """Generate the micro-structure

        Parameters
        ----------
        length : float
            length of the RVE
        width : float
            width of the RVE
        radius_mu : float
            mean radius of the fibers
        radius_std: float
            standard deviation of the fibers
        vol_req: float
            volume requirement of fiber material

        Returns
        -------
        volume_frac: float
            the actual volume fraction of the micro-structure.
        """

        microstructure_generator = CircleParticles(
            length=length,
            width=width,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
        )
        volume_frac = microstructure_generator.generate_rve()
        microstructure_generator.save_results()
        microstructure_generator.plot_rve(save_figure=True)

        return volume_frac

    def _save_data(self, responses) -> None:
        """save data to json file"""
        self.data["responses"] = responses
        working_folder = os.getcwd()
        with open("data.pickle", "wb") as file:
            pickle.dump(self.data, file)
        os.chdir(working_folder)

    def save_data(self, name: str = "data.pickle") -> None:
        working_folder = os.getcwd()
        with open(name, "wb") as file:
            pickle.dump(self.data, file)
        os.chdir(working_folder)
