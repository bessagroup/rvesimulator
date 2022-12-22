# import system packages
import os
import pickle

import pandas

import rvesimulator
from rvesimulator.microstructures.disks import CircleInclusion
from rvesimulator.simulators.abaqus_simulator import AbaqusSimulator
from rvesimulator.simulators.utils import create_dir


class SimulatorCaller:
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
            "post_script": "RVEPostProcess",
        }
        self.vol_req = None

    def run_simulation(self, data: dict = None) -> dict:
        """run simulation sequentially

        Parameters
        ----------
        data : dict, optional
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
            volume_frac = self.micro_structure_generation(
                length=self.rve_geometry["length"],
                width=self.rve_geometry["width"],
                radius=self.rve_geometry["radius"],
                vol_req=self.vol_req,
            )
            responses.at[ii, list(responses.keys())[0]] = volume_frac
            self._CompleteInformation()
            os.chdir(self.main_folder)
            abaqus_wrapper = AbaqusSimulator(
                sim_info=self.sim_info, folder_info=self.folder_info
            )
            abaqus_wrapper.Run()
            results = abaqus_wrapper.ReadBackResults()
            self.samples.at[ii, list(self.samples.keys())[2]] = results[
                list(self.samples.keys())[2]
            ]
            self.samples.at[ii, list(self.samples.keys())[3]] = results[
                list(self.samples.keys())[3]
            ]
            self._SaveData()

        return self.samples

    def update_sim_info(
        self,
        size: float = 0.048,
        radius: float = 0.003,
        vol_req: float = 0.30,
        mesh_partition: int = 30,
        loads: list = [0.05, 0.0, 0.0],
        simulation_time: float = 10.0,
        print_info: bool = False,
    ) -> None:
        """
        function to update the fixed simulation information
        Parameters
        ----------
        size (float): the size of the RVE, length=width=size
        radius (float): radius of the disks
        mesh_partition (int): the number of partition of every edge of the RVE
        loads(list): a list of loads [Exx, Eyy, Exy]
        simulation_time (float): total simulation time of Abaqus
        print_info(bool): a flag to print info or not

        Returns
        -------
        None
        """
        self.vol_req = vol_req
        self.rve_geometry = {"length": size, "width": size, "radius": radius}
        self.abaqus_paras = {
            "mesh_partition": mesh_partition,
            "loads": loads,
            "simulation_time": simulation_time,
        }
        if print_info:
            print(
                f"The general geometry information of RVE: {self.rve_geometry}"
            )
            print(f"The required volume fraction is: {self.volume_req}")
            print(
                f"The information of the Abaqus simulation : {self.abaqus_paras} \n"
            )

    def _update_sim_info(self, sample) -> None:
        """update the design variables"""
        self.sim_info.update(sample)

    @staticmethod
    def micro_structure_generation(
        length: float, width: float, radius: float, vol_req: float
    ) -> float:
        """Generate the micro-structure

        Parameters
        ----------
        length : float
            length of the RVE
        width : float
            width of the RVE
        radius : float
            radius of the RVE
        vol_req: float
            volume requirement of fiber material

        Returns
        -------
        volume_frac: float
            the actual volume fraction of the micro-structure.
        """

        microstructure_generator = CircleInclusion(
            length=length,
            width=width,
            radius=radius,
            vol_req=vol_req,
        )
        volume_frac = microstructure_generator.generate_rve()
        microstructure_generator.save_results()
        microstructure_generator.plot_rve()

        return volume_frac

    def _CompleteInformation(self) -> None:
        """
        This function is used to complete information for abaqus simulation

        Returns
        -------
        """

        file = "MicroStructureInfo.json"
        with open(file, "r") as f:
            location_info = json.load(f)

        self.sim_info = {
            "job_name": "ASCARVE",
            "location_information": location_info["location_information"],
            "Radius": location_info["Radius"],
            "LenStart": location_info["LenStart"],
            "LenEnd": location_info["LenEnd"],
            "WidStart": location_info["WidStart"],
            "WidEnd": location_info["WidEnd"],
            "mesh_partition": self.Abaqus_paras["mesh_partition"],
            "loads": self.Abaqus_paras["loads"],
            "simulation_time": self.Abaqus_paras["simulation_time"],
        }

    def _SaveData(self) -> None:
        """
        Function to save the simulation results into a Json file
        Returns
        -------

        """
        working_folder = os.getcwd()
        self.samples.to_json("doe.json", index=True)
        os.chdir(working_folder)
