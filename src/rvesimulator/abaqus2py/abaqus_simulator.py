"""
Module for calling abaqus simulation and post processing
"""
#                                                                       Modules
# =============================================================================
# Standard
import os
import pickle
import platform
import time
from pathlib import Path
from pickle import UnpicklingError
from typing import Dict, List

from ..additions.utils import (make_new_script, print_banner, remove_files,
                               write_json)
# local
from .simulator_interface import AssertInputs, Simulator

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class AbaqusSimulator(Simulator, AssertInputs):
    """Abaqus simulator

    Parameters
    ----------
    Simulator : class
        simulator interface
    AssertInputs: class
        assert inputs

    """

    def __init__(self,
                 sim_info: Dict,
                 folder_info: Dict) -> None:
        """initialization of abaqus simulator class

        Parameters
        ----------
        sim_info : Dict
            dict for defining abaqus simulation
        folder_info : Dict
            dict indicating abaqus scripts and working folders
        """

        self.sim_info = sim_info
        self.folder_info = folder_info
        # mandatory information
        self.is_job_name_in_sim_info(self.sim_info)
        self.is_script_path_in_folder_info(self.folder_info)
        self.is_mwd_in_folder_info(self.folder_info)
        self.job_name = sim_info["job_name"]
        self.main_dir = folder_info["main_dir"]
        self.script_path = folder_info["script_path"]

        # script information (not mandatory)
        self.current_dir = folder_info["current_dir"]
        self.sim_func = folder_info["sim_func"]
        self.sim_script = folder_info["sim_script"]
        self.post_func = folder_info["post_func"]
        self.post_script = folder_info["post_script"]
        # identify the platform of running abaqus
        self.platform = platform.system().lower()
        # define the name of json file
        self.abaqus_input_file = "sim_info.json"

    def pre_process(
        self,
        py_script: str,
        py_func: str,
    ) -> None:
        """execute the abaqus simulation

        Returns
        -------
        str
            simulation status
        """

        # hidden name information
        abaqus_py_script = str(self.job_name) + ".py"

        # write sim_info dict to a json file (every simulation has different
        # sim_info because some of them needs to be changed)
        write_json(sim_info=self.sim_info,
                   file_name=self.abaqus_input_file)

        # create new script for running abaqus simulation
        make_new_script(
            file_name=abaqus_py_script,
            script_path=self.script_path,
            script_name=py_script,
            func_name=py_func,
            input_file_name=self.abaqus_input_file,
        )

        # run abaqus simulation and submit the job
        print_banner("start .inp generation ")

        # call abaqus no-gui via system command
        self._call_abaqus_no_gui(file_name=abaqus_py_script)

    def submit_job(self, num_cpu: int = 4) -> None:

        print_banner("submit .inp to abaqus solver")

        file_name = self.job_name + ".inp"
        if "subroutine_path" not in self.sim_info.keys():
            command = f"abaqus job={file_name} cpus={num_cpu}  -interactive"
        else:
            command = f"abaqus job={file_name} cpus={num_cpu} user={self.sim_info['subroutine_path']} -interactive"
        start_time = time.perf_counter()
        os.system(command)
        end_time = time.perf_counter()
        print(f"abaqus solver finished with  :{(end_time - start_time):2f} s")

    def post_process(self,
                     post_py_script: str,
                     post_py_func: str,
                     delete_odb: bool = False) -> None:
        """post process to get the preliminary results

        Parameters
        ----------
        delete_odb : bool, optional
            delete odb file to save memory, by default True
        """
        print_banner("start post processing ")
        if post_py_func is not None and post_py_script is not None:
            # path with the python-script
            post_process_script = "getResults.py"

            # make new script for post processing
            make_new_script(
                file_name=post_process_script,
                script_path=self.script_path,
                script_name=post_py_script,
                func_name=post_py_func,
                input_file_name=self.abaqus_input_file,
            )

            # call abaqus no-gui via system command
            self._call_abaqus_no_gui(file_name=post_process_script)

            # remove files that influence the simulation process
            remove_files(directory=os.getcwd())

        # remove the odb file to save memory
        if delete_odb:
            remove_files(directory=os.getcwd(), file_types=[".odb"])

    def read_back_results(self, file_name: str = "results.pkl") -> Dict:
        """read back the results from the generated pickle file

        Parameters
        ----------
        file_name : str, optional
            file name, by default "results.pkl"

        Returns
        -------
        dict
            a dict containing the results
        """

        try:
            with open(file_name, "rb") as fd:
                results = pickle.load(fd, fix_imports=True, encoding="latin1")
        except UnpicklingError:
            # fix issue of windows system, if the simulation is ran on windows
            # then the pickle file is not been read via above method. In this
            # case, the following code is used to work around this issue
            # TODO: find a better way to fix this issue
            content = ''
            outsize = 0
            with open(file_name, 'rb') as infile:
                content = infile.read()
            with open(file_name, 'wb') as output:
                for line in content.splitlines():
                    outsize += len(line) + 1
                    output.write(line + str.encode('\n'))
            # open the file again
            with open(file_name, "rb") as fd:
                results = pickle.load(fd, fix_imports=True, encoding="latin1")
        return results

    def _call_abaqus_no_gui(
        self,
        file_name: str,
    ) -> None:
        """run abaqus simulation

        Returns
        -------
        str
            simulation status

        Raises
        ------
        NotImplementedError
            platform not be implemented
        """

        # count time for simulation
        start_time = time.perf_counter()
        command = "abaqus cae noGUI=" + str(file_name) + " -mesa"
        os.system(command)
        end_time = time.perf_counter()

        # detect the existence of the inp file
        inp_exist = Path(self.job_name + ".inp").exists()
        odb_exist = Path(self.job_name + ".odb").exists()

        if inp_exist and odb_exist:
            print(
                f"post process finished with  :{(end_time - start_time):2f} s")
        elif not odb_exist and inp_exist:
            print(
                f"inp generation finished with  :{(end_time - start_time):2f} s")
