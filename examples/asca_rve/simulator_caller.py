# import system packages
import json
import os
import pickle

# local functions
import rvesimulator
from rvesimulator.microstructures.homo_radius_disks import CircleInclusion
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
        self.update_sim_info()

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
            # update the geometry info for microstructure
            self._update_sample_info(sample=samples[ii])

            # generating the micro-structure
            vol_frac = self.micro_structure_generation(
                length=self.rve_geometry["length"],
                width=self.rve_geometry["width"],
                radius=self.rve_geometry["radius"],
                vol_req=self.vol_req,
            )
            list_iter = list(responses.keys())
            if "vol_frac" in list_iter:
                responses.at[ii, "vol_frac"] = vol_frac
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
        radius: float = 0.003,
        vol_req: float = 0.30,
        mesh_partition: int = 30,
        loads: list = [0.05, 0.0, 0.0],
        simulation_time: float = 10.0,
        print_info: bool = False,
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
        self.rve_geometry = {"length": size, "width": size, "radius": radius}
        self.abaqus_paras = {
            "mesh_partition": mesh_partition,
            "loads": loads,
            "simulation_time": simulation_time,
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
            "radius": location_info["radius"],
            "len_start": location_info["len_start"],
            "len_end": location_info["len_end"],
            "wid_start": location_info["wid_start"],
            "wid_end": location_info["wid_end"],
            "mesh_partition": self.abaqus_paras["mesh_partition"],
            "loads": self.abaqus_paras["loads"],
            "simulation_time": self.abaqus_paras["simulation_time"],
        }

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
            second_heuristic=False,
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
