# import system packages
import os
import pickle
import subprocess
import time

# import local functions
from rvesimulator.simulators.utils import (
    create_dir,
    kill_abaqus_process,
    make_new_script,
    print_banner,
    remove_files,
    write_json,
)

from rvesimulator.simulators.simulator_interface import Simulator


class AbaqusSimulator(Simulator):
    """
    This class is used for running Abaqus simulation and getting responses
    """

    def __init__(self, sim_info: dict, folder_info: dict) -> None:
        """intialization of abaqus simulator

        Parameters
        ----------
        sim_info : dict
            simulation information for abaqus
        folder_info : dict
            folders information
        """

        self.sim_info = sim_info
        self.folder_info = folder_info
        self.job_name = sim_info["job_name"]
        self.main_work_directory = folder_info["main_work_directory"]
        self.current_work_directory = folder_info["current_work_directory"]
        self.script_path = folder_info["script_path"]
        self.sim_path = folder_info["sim_path"]
        self.sim_script = folder_info["sim_script"]
        self.post_path = folder_info["post_path"]
        self.post_script = folder_info["post_script"]

        # define the name of json file
        self.sim_info_name = "sim_info.json"

    def execute(self) -> None:
        """execute simulation"""

        # hidden name information
        new_python_filename = "abqScript.py"

        # folder operations
        new_path = create_dir(
            current_folder=self.main_work_directory,
            dirname=self.current_work_directory,
        )

        # change work directory
        os.chdir(new_path)
        print("Current working directory: {0}".format(os.getcwd()))
        # output the sim_info dict to a json file
        write_json(sim_info=self.sim_info, filename=self.sim_info_name)
        #
        make_new_script(
            new_file_name=new_python_filename,
            folder_info=self.folder_info,
            status="simulation",
            sim_info_name=self.sim_info_name,
        )

        # begin to run abaqus simulation and submit the job to get the .odb file
        print_banner("START ABAQUS ANALYSIS")
        # the first step is write the design variables into a .txt file
        start_time = time.time()
        command = "abaqus cae noGUI=" + str(new_python_filename) + " -mesa"
        self._run_abaqus_simulation(command)
        self._remove_files()
        time_end = time.time()
        print(f"time cost of this iteraion: {time_end - start_time}")
        # change the work directory back the main one

    def post_process(self, delete_odb: bool = True) -> None:
        """post process

        Parameters
        ----------
        delete_odb : bool, optional
            delete the odb file or not, by default True
        """
        print_banner("START ABAQUS POST ANALYSIS")
        # path with the python-script
        new_python_filename = "get_results.py"
        make_new_script(
            new_file_name=new_python_filename,
            folder_info=self.folder_info,
            status="post_process",
            job_name=self.job_name,
        )
        command = "abaqus cae noGUI=" + str(new_python_filename) + " -mesa"
        os.system(command)
        self._remove_files()

        if delete_odb:
            odb_file = self.job_name + ".odb"
            if os.path.exists(odb_file):
                os.remove(odb_file)

    def read_back_results(self) -> dict:
        """read the simulation results back to python

        Returns
        -------
        results : dict
            simulation results
        """

        with open("results.p", "rb") as fd:
            results = pickle.load(fd, fix_imports=True, encoding="latin1")
        os.chdir(self.main_work_directory)

        return results

    def _run_abaqus_simulation(self, command: str) -> None:
        """This function is used to run abaqus simulation

        Parameters
        ----------
        command : str
            system command to run abaqus
        """

        my_cmd_job = command
        proc = subprocess.Popen(my_cmd_job, shell=True)
        start_time = time.time()
        time.sleep(20.0)
        while True:
            # in this part: check the job is finish or not !
            time.sleep(20.0 - ((time.time() - start_time) % 20.0))
            msg_name = self.job_name + ".msg"
            file = open(msg_name)
            word1 = "THE ANALYSIS HAS BEEN COMPLETED"
            if word1 in file.read():
                print("Simulation successfully finished! \n")
                proc.kill()
                kill_abaqus_process()
                break
            end_time = time.time()
            print(f"the simulation time is :{end_time - start_time} !")

    def _remove_files(self) -> None:
        """remove files to save memory"""
        log_file = self.job_name + ".log"
        if os.path.exists(log_file):
            os.remove(log_file)
        lck_file = self.job_name + ".lck"
        if os.path.exists(lck_file):
            os.remove(lck_file)
        directory = os.getcwd()

        list_of_removed_filetypes = [
            ".SMABulk",
            ".rec",
            ".SMAFocus",
            ".exception",
            ".simlog",
        ]
        for filetype in list_of_removed_filetypes:
            remove_files(directory, filetype)
