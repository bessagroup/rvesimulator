# import system packages
import json
import os
import pickle

import numpy as np

# local functions
import rvesimulator
from rvesimulator.microstructures.heter_radius_circles import (
    HeterCircleInclusion,
)
from rvesimulator.simulators.abaqus_simulator import AbaqusSimulator
from rvesimulator.simulators.utils import create_dir


class CooperativeRVE:
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
            "sim_path": "scriptbase.pnas_rve",
            "sim_script": "PnasCompositeRVE",
            "post_path": "scriptbase.postprocess",
            "post_script": "RVEPostProcess2D",
        }
        self.vol_req = None
        # self.update_sim_info(yield_stress=0.5, a=0.5)

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
                # remove the result
                self.remove_file(file_name="results.p")
                # remove the micro-structure file
                if self.update_micro_structure is True:
                    self.remove_file(file_name="micro_structure_info.json")
            else:
                self.folder_info["current_work_directory"] = "data"
                new_path = create_dir(
                    current_folder=self.folder_info["main_work_directory"],
                    dirname=self.folder_info["current_work_directory"],
                )
                os.chdir(new_path)
                # remove the result
                self.remove_file(file_name="results.p")
                # remove the micro-structure file
                if self.update_micro_structure is True:
                    self.remove_file(file_name="micro_structure_info.json")

            # update the geometry info for microstructure
            self._update_sample_info(sample=samples[ii])
            if os.path.isfile("micro_structure_info.json"):
                print("micro-structure file already exist\n")
            else:
                self.vol_frac = self.micro_structure_generation(
                    length=self.rve_geometry["length"],
                    width=self.rve_geometry["width"],
                    radius_mu=self.rve_geometry["radius_mu"],
                    radius_std=self.rve_geometry["radius_std"],
                    vol_req=self.vol_req,
                    seed=self.rve_geometry["seed"],
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
        radius_std: float = 0,
        vol_req: float = 0.30,
        mesh_partition: int = 100,
        loads: list = [0.02, 0.02, 0.02],
        loads_path: list = None,
        E_matrix: float = 100.0,
        Pr_matrix: float = 0.30,
        hardening_law: str = "linear",
        E_fiber: float = 1.0,
        Pr_fiber: float = 0.19,
        time_period: float = 1.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        update_micro_structure: bool = False,
        print_info: bool = False,
        seed: int = None,
        **kwarg,
    ) -> None:
        """update some default information

        Parameters
        ----------
        size : float, optional
            size of RVE, by default 0.048
        radius : float, optional
            radius of particles, by default 0.003
        vol_req : float, optional
            required volume fraction, by default 0.30
        mesh_partition : int, optional
            mesh partition, by default 30
        loads : list, optional
            boundary of loads, by default [0.05, 0.0, 0.0]
        loads_path : list, optional
            loads path, by default None
        E_matrix : float, optional
            Young's modulus of matrix material, by default 100.0
        Pr_matrix : float, optional
            poission ratio of matrix material, by default 0.30
        matrix_yield_law : str, optional
            plasticity law of matrix material, by default "Von_mises"
        E_fiber : float, optional
            Young's modulus of fiber material, by default 1.0
        Pr_fiber : float, optional
            poission ratio of fiber material, by default 0.19
        time_period : float, optional
            simulation time, by default 1.0
        print_info : bool, optional
            a flag to indicate print the default/update simulation information
            , by default False

        Raises
        ------
        KeyError
            key error for matrix material plasticity law definition string
        """
        # to know if I want to updatet the micro-structure file or not
        self.update_micro_structure = update_micro_structure
        # print the info to screen
        if update_micro_structure is True:
            print("Micro-structure file will be updated \n")
        else:
            print("Micro-structure file will not be updated \n")
        # generate the yield function table for matrix material
        yield_table_matrix = self.hardening_law(
            law_type=hardening_law, **kwarg
        )
        # define the parameters for micro-structure
        self.vol_req = vol_req
        self.rve_geometry = {
            "length": size,
            "width": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "seed": seed,
        }
        self.abaqus_paras = {
            "mesh_partition": mesh_partition,
            "loads": loads,
            "time_period": time_period,
            "loads_path": loads_path,
            "E_matrix": E_matrix,
            "Pr_matrix": Pr_matrix,
            "yield_table_matrix": yield_table_matrix,
            "E_fiber": E_fiber,
            "Pr_fiber": Pr_fiber,
            "num_cpu": num_cpu,
            "platform": platform,
        }
        if print_info:
            print(f"geometry information of RVE: {self.rve_geometry}")
            print(f"vol_req is: {self.vol_req}")
            print(f"Info of Abaqus simulation : {self.abaqus_paras} \n")

    def _update_sample_info(self, sample: dict) -> None:
        """_summary_

        Parameters
        ----------
        sample : dict
            samples at current iteration
        """

        # update 'vol_req'
        if "vol_req" in sample.keys():
            self.vol_req = sample["vol_req"]

        # update 'size'
        if "size" in sample.keys():
            self.rve_geometry["length"] = sample["size"]
            self.rve_geometry["width"] = sample["size"]

        # update 'radius'
        if "radius_mu" in sample.keys():
            self.rve_geometry["radius_mu"] = sample["radius_mu"]

        if "radius_std" in sample.keys():
            self.rve_geometry["radius_std"] = sample["radius_std"]

        # update 'mesh'
        if "mesh_partition" in sample.keys():
            self.abaqus_paras["mesh_partition"] = sample["mesh_partition"]
        # update 'loads path'
        if "loads_path" in sample.keys():
            self.abaqus_paras["loads_path"] = sample["loads_path"]

        if "E_matrix" in sample.keys():
            self.abaqus_paras["E_matrix"] = sample["E_matrix"]

        if "Pr_matrix" in sample.keys():
            self.abaqus_paras["Pr_matrix"] = sample["Pr_matrix"]

        if "E_fiber" in sample.keys():
            self.abaqus_paras["E_fiber"] = sample["E_fiber"]

        if "Pr_fiber" in sample.keys():
            self.abaqus_paras["Pr_fiber"] = sample["Pr_fiber"]

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
            "job_name": "pnas_composite",
            "location_information": location_info["location_information"],
            "radius_mu": location_info["radius_mu"],
            "radius_std": location_info["radius_std"],
            "len_start": location_info["len_start"],
            "len_end": location_info["len_end"],
            "wid_start": location_info["wid_start"],
            "wid_end": location_info["wid_end"],
            "mesh_partition": self.abaqus_paras["mesh_partition"],
            "loads": self.abaqus_paras["loads"],
            "loads_path": self.abaqus_paras["loads_path"],
            "E_matrix": self.abaqus_paras["E_matrix"],
            "Pr_matrix": self.abaqus_paras["Pr_matrix"],
            "yield_table_matrix": self.abaqus_paras["yield_table_matrix"],
            "E_fiber": self.abaqus_paras["E_fiber"],
            "Pr_fiber": self.abaqus_paras["Pr_fiber"],
            "time_period": self.abaqus_paras["time_period"],
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
        seed: int = None,
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

        microstructure_generator = HeterCircleInclusion(
            length=length,
            width=width,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )
        volume_frac = microstructure_generator.generate_rve()
        microstructure_generator.save_results()
        microstructure_generator.plot_rve(save_figure=False)

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

    @staticmethod
    def hardening_law(law_type: str = "linear", **kwarg) -> list:
        """yield crterion $\sigma_y = factor_1 + factor_2 \times exp(\epsilon)^factor_3$

        Parameters
        ----------
        law_type: str
            a str that indicates the hardening law

        Returns
        -------
        yield_table : list
            a list contains the yield table
        """
        yield_table = np.zeros((101, 2))
        yield_table[:, 1] = np.linspace(0, 1, 101)
        if law_type == "linear":
            yield_stress = kwarg["yield_stress"]
            a = kwarg["a"]
            yield_table[:, 0] = yield_stress + a * yield_table[:, 1]
            yield_table[-1, 1] = 10.0
            yield_table[-1, 0] = yield_stress + a * yield_table[-1, 1]

        elif law_type == "swift":
            yield_stress = kwarg["yield_stress"]
            a = kwarg["a"]
            b = kwarg["b"]
            yield_table[:, 0] = yield_stress + a * (yield_table[:, 1]) ** b
            yield_table[-1, 1] = 10.0
            yield_table[-1, 0] = yield_stress + a * (yield_table[-1, 1]) ** b

        elif law_type == "ramberg":
            yield_stress = kwarg["yield_stress"]
            a = kwarg["a"]
            b = kwarg["b"]
            yield_table[:, 0] = yield_stress * (
                1 + a * (yield_table[:, 1])
            ) ** (1 / b)
            yield_table[-1, 1] = 10.0
            yield_table[-1, 0] = yield_stress * (
                1 + a * (yield_table[-1, 1])
            ) ** (1 / b)
        else:
            raise KeyError("This hardening law is not defined \n")

        yield_table = yield_table.T
        return yield_table.tolist()

    @staticmethod
    def remove_file(file_name: str = "micro_structure_info.json") -> None:
        """remove file

        Parameters
        ----------
        file_name : str, optional
            name of the file, by default "micro_structure_info.json"
        """

        if os.path.exists(file_name):
            print(f"remove {file_name} successfully\n")
            os.remove(file_name)
        else:
            print(f"{file_name} do not exist\n")


class CooperativeRVESwap:
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
            "sim_path": "scriptbase.pnas_rve",
            "sim_script": "PnasCompositeRVESwap",
            "post_path": "scriptbase.postprocess",
            "post_script": "RVEPostProcess2D",
        }
        self.vol_req = None
        # self.update_sim_info(yield_stress=0.5, a=0.5)

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
                # remove the result
                self.remove_file(file_name="results.p")
                # remove the micro-structure file
                if self.update_micro_structure is True:
                    self.remove_file(file_name="micro_structure_info.json")
            else:
                self.folder_info["current_work_directory"] = "data"
                new_path = create_dir(
                    current_folder=self.folder_info["main_work_directory"],
                    dirname=self.folder_info["current_work_directory"],
                )
                os.chdir(new_path)
                # remove the result
                self.remove_file(file_name="results.p")
                # remove the micro-structure file
                if self.update_micro_structure is True:
                    self.remove_file(file_name="micro_structure_info.json")

            # update the geometry info for microstructure
            self._update_sample_info(sample=samples[ii])
            if os.path.isfile("micro_structure_info.json"):
                print("micro-structure file already exist\n")
            else:
                self.vol_frac = self.micro_structure_generation(
                    length=self.rve_geometry["length"],
                    width=self.rve_geometry["width"],
                    radius_mu=self.rve_geometry["radius_mu"],
                    radius_std=self.rve_geometry["radius_std"],
                    vol_req=self.vol_req,
                    seed=self.rve_geometry["seed"],
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
        radius_std: float = 0,
        vol_req: float = 0.30,
        mesh_partition: int = 100,
        loads: list = [0.02, 0.02, 0.02],
        loads_path: list = None,
        E_matrix: float = 1.0,
        Pr_matrix: float = 0.19,
        hardening_law: str = "linear",
        E_fiber: float = 100.0,
        Pr_fiber: float = 0.30,
        time_period: float = 1.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        update_micro_structure: bool = False,
        print_info: bool = False,
        seed: int = None,
        **kwarg,
    ) -> None:
        """update some default information

        Parameters
        ----------
        size : float, optional
            size of RVE, by default 0.048
        radius : float, optional
            radius of particles, by default 0.003
        vol_req : float, optional
            required volume fraction, by default 0.30
        mesh_partition : int, optional
            mesh partition, by default 30
        loads : list, optional
            boundary of loads, by default [0.05, 0.0, 0.0]
        loads_path : list, optional
            loads path, by default None
        E_matrix : float, optional
            Young's modulus of matrix material, by default 100.0
        Pr_matrix : float, optional
            poission ratio of matrix material, by default 0.30
        matrix_yield_law : str, optional
            plasticity law of matrix material, by default "Von_mises"
        E_fiber : float, optional
            Young's modulus of fiber material, by default 1.0
        Pr_fiber : float, optional
            poission ratio of fiber material, by default 0.19
        time_period : float, optional
            simulation time, by default 1.0
        print_info : bool, optional
            a flag to indicate print the default/update simulation information
            , by default False

        Raises
        ------
        KeyError
            key error for matrix material plasticity law definition string
        """
        # to know if I want to updatet the micro-structure file or not
        self.update_micro_structure = update_micro_structure
        # print the info to screen
        if update_micro_structure is True:
            print("Micro-structure file will be updated \n")
        else:
            print("Micro-structure file will not be updated \n")
        # generate the yield function table for matrix material
        yield_table_fiber = self.hardening_law(law_type=hardening_law, **kwarg)
        # define the parameters for micro-structure
        self.vol_req = vol_req
        self.rve_geometry = {
            "length": size,
            "width": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "seed": seed,
        }
        self.abaqus_paras = {
            "mesh_partition": mesh_partition,
            "loads": loads,
            "time_period": time_period,
            "loads_path": loads_path,
            "E_matrix": E_matrix,
            "Pr_matrix": Pr_matrix,
            "yield_table_fiber": yield_table_fiber,
            "E_fiber": E_fiber,
            "Pr_fiber": Pr_fiber,
            "num_cpu": num_cpu,
            "platform": platform,
        }
        if print_info:
            print(f"geometry information of RVE: {self.rve_geometry}")
            print(f"vol_req is: {self.vol_req}")
            print(f"Info of Abaqus simulation : {self.abaqus_paras} \n")

    def _update_sample_info(self, sample: dict) -> None:
        """_summary_

        Parameters
        ----------
        sample : dict
            samples at current iteration
        """

        # update 'vol_req'
        if "vol_req" in sample.keys():
            self.vol_req = sample["vol_req"]

        # update 'size'
        if "size" in sample.keys():
            self.rve_geometry["length"] = sample["size"]
            self.rve_geometry["width"] = sample["size"]

        # update 'radius'
        if "radius_mu" in sample.keys():
            self.rve_geometry["radius_mu"] = sample["radius_mu"]

        if "radius_std" in sample.keys():
            self.rve_geometry["radius_std"] = sample["radius_std"]

        # update 'mesh'
        if "mesh_partition" in sample.keys():
            self.abaqus_paras["mesh_partition"] = sample["mesh_partition"]
        # update 'loads path'
        if "loads_path" in sample.keys():
            self.abaqus_paras["loads_path"] = sample["loads_path"]

        if "E_matrix" in sample.keys():
            self.abaqus_paras["E_matrix"] = sample["E_matrix"]

        if "Pr_matrix" in sample.keys():
            self.abaqus_paras["Pr_matrix"] = sample["Pr_matrix"]

        if "E_fiber" in sample.keys():
            self.abaqus_paras["E_fiber"] = sample["E_fiber"]

        if "Pr_fiber" in sample.keys():
            self.abaqus_paras["Pr_fiber"] = sample["Pr_fiber"]

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
            "job_name": "pnas_composite",
            "location_information": location_info["location_information"],
            "radius_mu": location_info["radius_mu"],
            "radius_std": location_info["radius_std"],
            "len_start": location_info["len_start"],
            "len_end": location_info["len_end"],
            "wid_start": location_info["wid_start"],
            "wid_end": location_info["wid_end"],
            "mesh_partition": self.abaqus_paras["mesh_partition"],
            "loads": self.abaqus_paras["loads"],
            "loads_path": self.abaqus_paras["loads_path"],
            "E_matrix": self.abaqus_paras["E_matrix"],
            "Pr_matrix": self.abaqus_paras["Pr_matrix"],
            "yield_table_fiber": self.abaqus_paras["yield_table_fiber"],
            "E_fiber": self.abaqus_paras["E_fiber"],
            "Pr_fiber": self.abaqus_paras["Pr_fiber"],
            "time_period": self.abaqus_paras["time_period"],
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
        seed: int = None,
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

        microstructure_generator = HeterCircleInclusion(
            length=length,
            width=width,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )
        volume_frac = microstructure_generator.generate_rve()
        microstructure_generator.save_results()
        # microstructure_generator.plot_rve(save_figure=False)

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

    @staticmethod
    def hardening_law(law_type: str = "linear", **kwarg) -> list:
        """yield crterion $\sigma_y = factor_1 + factor_2 \times exp(\epsilon)^factor_3$

        Parameters
        ----------
        law_type: str
            a str that indicates the hardening law

        Returns
        -------
        yield_table : list
            a list contains the yield table
        """
        yield_table = np.zeros((101, 2))
        yield_table[:, 1] = np.linspace(0, 1, 101)
        if law_type == "linear":
            yield_stress = kwarg["yield_stress"]
            a = kwarg["a"]
            yield_table[:, 0] = yield_stress + a * yield_table[:, 1]
            yield_table[-1, 1] = 10.0
            yield_table[-1, 0] = yield_stress + a * yield_table[-1, 1]

        elif law_type == "swift":
            yield_stress = kwarg["yield_stress"]
            a = kwarg["a"]
            b = kwarg["b"]
            yield_table[:, 0] = yield_stress + a * (yield_table[:, 1]) ** b
            yield_table[-1, 1] = 10.0
            yield_table[-1, 0] = yield_stress + a * (yield_table[-1, 1]) ** b

        elif law_type == "ramberg":
            yield_stress = kwarg["yield_stress"]
            a = kwarg["a"]
            b = kwarg["b"]
            yield_table[:, 0] = yield_stress * (
                1 + a * (yield_table[:, 1])
            ) ** (1 / b)
            yield_table[-1, 1] = 10.0
            yield_table[-1, 0] = yield_stress * (
                1 + a * (yield_table[-1, 1])
            ) ** (1 / b)
        else:
            raise KeyError("This hardening law is not defined \n")

        yield_table = yield_table.T
        return yield_table.tolist()

    @staticmethod
    def remove_file(file_name: str = "micro_structure_info.json") -> None:
        """remove file

        Parameters
        ----------
        file_name : str, optional
            name of the file, by default "micro_structure_info.json"
        """

        if os.path.exists(file_name):
            print(f"remove {file_name} successfully\n")
            os.remove(file_name)
        else:
            print(f"{file_name} do not exist\n")
