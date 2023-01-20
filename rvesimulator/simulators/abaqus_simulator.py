#                                                                       Modules
# =============================================================================
# Standard
import json
import os
import pickle
import subprocess
import time

# Local
from rvesimulator.simulators.simulator_interface import Simulator
from rvesimulator.simulators.utils import create_dir, write_json

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class AbaqusSimulator(Simulator):
    """
    running Abaqus simulation and getting responses
    """

    def __init__(self, sim_info: dict, folder_info: dict) -> None:

        """intialization of abaqus simulator

        Parameters
        ----------
        sim_info : dict
            simulation information for abaqus
        folder_info : dict
            folders information
        plaform : bool
            platform of running simulation
        """

        self.sim_info = sim_info
        self.folder_info = folder_info

        # details of variables
        self.job_name = sim_info["job_name"]
        self.main_work_directory = folder_info["main_work_directory"]
        self.current_work_directory = folder_info["current_work_directory"]

        # script information
        self.script_path = folder_info["script_path"]
        self.sim_path = folder_info["sim_path"]
        self.sim_script = folder_info["sim_script"]

        # only used on 'ubuntu' platform
        self.post_path = folder_info["post_path"]
        self.post_script = folder_info["post_script"]
        self.platform = sim_info["platform"]

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
        # create new script for running abaqus simulation
        self.make_new_script(
            new_file_name=new_python_filename,
            folder_info=self.folder_info,
            status="simulation",
            sim_info_name=self.sim_info_name,
        )

        # run abaqus simulation and submit the job
        self.print_banner("START ABAQUS ANALYSIS")
        start_time = time.time()
        command = "abaqus cae noGUI=" + str(new_python_filename) + " -mesa"
        if self.platform == "cluster":
            os.system(command)
        elif self.platform == "ubuntu":
            self._run_abaqus_simulation(command)
            self._remove_files()
        elif self.platform == "windows":
            raise NotImplementedError("it is not implement yet \n")
        else:
            raise NotImplementedError("it is not implement yet \n")
        time_end = time.time()
        print(f"time cost of this iteraion: {time_end - start_time}")

    def post_process(self, delete_odb: bool = True) -> None:
        """post process

        Parameters
        ----------
        delete_odb : bool, optional
            delete the odb file or not, by default True
        """
        if self.platform == "ubuntu":
            self.print_banner("START ABAQUS POST ANALYSIS")
            # path with the python-script
            new_python_filename = "get_results.py"
            self.make_new_script(
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
        elif self.platform == "cluster":
            pass

        elif self.platform == "windows":
            pass

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
        """This function is used to run abaqus simulation, only used when
        running simulation on the "ubuntu" platform

        Parameters
        ----------
        command : str
            system command to run abaqus
        """

        my_cmd_job = command
        proc = subprocess.Popen(my_cmd_job, shell=True)
        start_time = time.time()
        time.sleep(30.0)
        while True:
            # in this part: check the job is finish or not !
            time.sleep(20.0 - ((time.time() - start_time) % 20.0))
            msg_name = self.job_name + ".msg"
            file = open(msg_name)
            word1 = "THE ANALYSIS HAS BEEN COMPLETED"
            if word1 in file.read():
                print("Simulation successfully finished! \n")
                proc.kill()
                self.kill_abaqus_process()
                break
            end_time = time.time()
            print(f"simulation time :{(end_time - start_time):2f} s")

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
            self.remove_files(directory, filetype)

    @staticmethod
    def make_new_script(
        new_file_name: str,
        folder_info: dict,
        status: str = "simulation",
        sim_info_name: str = "sim_info.json",
        job_name: str = "Job-1",
    ) -> None:
        """make a new python script to bridge Python and abaqus

        Parameters
        ----------
        new_file_name : str
            new file name
        folder_info : dict
            folder where the new script will be created
        status : str, optional
            status of simulation, could be "simulation" and "post_process",
            by default "simulation"
        sim_info_name : str, optional
            the name of file contains input variables for abaqus simulation,
            by default "sim_info.json"
        job_name : str, optional
            job name of simulation,  by default "Job-1"
        """

        if status == "simulation":
            with open(new_file_name, "w") as file:
                file.write("import os \n")
                file.write("import sys \n")
                file.write("import json \n")
                file.write(
                    "sys.path.extend(['"
                    + str(folder_info["script_path"])
                    + "']) \n"
                )
                file.write(
                    "from "
                    + str(folder_info["sim_path"])
                    + " import "
                    + str(folder_info["sim_script"])
                    + "\n"
                )
                line = "file = '" + str(sim_info_name) + "' \n"
                file.write(line)
                file.write("with open(file, 'r') as f:\n")
                file.write("	dict = json.load(f)\n")
                file.write(str(folder_info["sim_script"]) + "(dict)\n")
            file.close()
        elif status == "post_process":
            with open(new_file_name, "w") as file:
                file.write("import os\n")
                file.write("import sys\n")
                file.write(
                    "sys.path.extend(['"
                    + str(folder_info["script_path"])
                    + "']) \n"
                )
                file.write(
                    "from "
                    + str(folder_info["post_path"])
                    + " import "
                    + str(folder_info["post_script"])
                    + "\n"
                )
                file.write(
                    str(folder_info["post_script"])
                    + "('"
                    + str(job_name)
                    + "')\n"
                )
            file.close()
        else:
            raise KeyError("the process is not needed for a new script \n")

    @staticmethod
    def kill_abaqus_process() -> None:
        """kill the simulation process because it can not stop properly"""
        name_1 = "pkill standard"
        aa = os.system(name_1)
        name_2 = "pkill ABQcaeK"
        bb = os.system(name_2)
        name_3 = "pkill SMAPython"
        cc = os.system(name_3)
        print(aa + bb + cc)

    @staticmethod
    def print_banner(message: str, sign="#", length=50) -> None:
        """print banner

        Parameters
        ----------
        message : str
            string output to the screen
        sign : str, optional
            pattern, by default "#"
        length : int, optional
            length of output, by default 50
        """
        print(sign * length)
        print(
            sign * ((length - len(message) - 2) // 2)
            + " "
            + message
            + " "
            + sign * ((length - len(message) - 2) // 2)
        )
        print(sign * length)

    @staticmethod
    def remove_files(directory: str, file_type: str) -> None:
        """remove files

        Parameters
        ----------
        directory : _type_
            working directory
        file_type : str
            file type
        """
        files_in_directory = os.listdir(directory)
        filtered_files = [
            file for file in files_in_directory if file.endswith(file_type)
        ]
        for file in filtered_files:
            path_to_file = os.path.join(directory, file)
            if os.path.exists(path_to_file):
                os.remove(path_to_file)
