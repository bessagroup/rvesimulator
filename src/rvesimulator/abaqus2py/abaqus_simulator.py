"""
Module for calling abaqus simulation and post processing
"""
#                                                                       Modules
# =============================================================================
# Standard
import os
import pickle
import subprocess
import time
from math import inf
from pathlib import Path
from pickle import UnpicklingError

# local
from .simulator_interface import AssertInputs, Simulator
from .utils import print_banner, write_json

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

    Examples
    --------
    >>> from rvesimulator.abaqus2py.simulator_interface import AbaqusSimulator
    >>> sim_info = {}
    >>> folder_info = {}
    >>> simulator = AbaqusSimulator(sim_info=sim_info, folder_info=folder_info)
    >>> simulator.run()
    >>> simulator.read_back_results()
    """

    def __init__(self, sim_info: dict, folder_info: dict) -> None:
        """initialization of abaqus simulator class

        Parameters
        ----------
        sim_info : dict
            dict for defining abaqus simulation
        folder_info : dict
            dict indicating abaqus scripts and working folders
        """

        self.sim_info = sim_info
        self.folder_info = folder_info
        # to assert the needed information is proper defined
        self.is_inputs_proper_defined(
            folder_info=folder_info, sim_info=sim_info
        )

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
        #
        self.platform = sim_info["platform"]

        # define the name of json file
        self.abaqus_input_file = "sim_info.json"

    def execute(
        self,
        max_time: float = inf,
        sleep_time: float = 20.0,
        refresh_time: float = 5.0,
    ) -> str:
        """execute the abaqus simulation

        Parameters
        ----------
        max_time : float, optional
            maximum simulation time, kill the simulation when time reaches the
            given maximum simulation time, by default inf
        sleep_time : float, optional
            sleep time, let abaqus run for 20s first, by default 20.0
        refresh_time : float, optional
            detect if the simulation is finished or not every given refresh
            time, by default 5.0

        Returns
        -------
        str
            simulation status
        """

        # hidden name information
        abaqus_py_script = "abaScript.py"

        # write sim_info dict to a json file
        write_json(sim_info=self.sim_info, file_name=self.abaqus_input_file)

        # create new script for running abaqus simulation
        self.make_new_script(
            file_name=abaqus_py_script,
            status="simulation",
        )

        # run abaqus simulation and submit the job
        print_banner("start abaqus simulation")
        # system command for running abaqus
        command = "abaqus cae noGUI=" + str(abaqus_py_script) + " -mesa"

        # run (flag is used to identify if the simulation is finished or killed
        # because of overtime)
        self.flag = self._run_abaqus_simulation(
            command=command,
            max_time=max_time,
            sleep_time=sleep_time,
            refresh_time=refresh_time,
        )

        return self.flag

    def post_process(self, delete_odb: bool = False) -> None:
        """post process to get the preliminary results

        Parameters
        ----------
        delete_odb : bool, optional
            delete odb file to save memory, by default True
        """

        if self.platform == "ubuntu":
            if self.flag == "error":
                self.job_name = "error"
            print_banner("abaqus post analysis")
            # path with the python-script
            post_process_script = "getResults.py"
            self.make_new_script(
                file_name=post_process_script,
                status="post_process",
            )
            command = "abaqus cae noGUI=" + str(post_process_script) + " -mesa"
            os.system(command)

            # remove files that influence the simulation process
            self.remove_files(directory=os.getcwd())

        # remove the odb file to save memory
        if delete_odb:
            self.remove_files(directory=os.getcwd(), file_types=[".odb"])

    def read_back_results(self, file_name: str = "results.p") -> dict:
        """read back the results from the generated pickle file

        Parameters
        ----------
        file_name : str, optional
            file name, by default "results.p"

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

    def _run_abaqus_simulation(
        self,
        command: str,
        max_time: float = inf,
        sleep_time: float = 20.0,
        refresh_time: float = 5.0,
    ) -> str:
        """run abaqus simulation

        Parameters
        ----------
        command : str
            system command for running abaqus
        max_time : float, optional
            maximum simulation time, by default inf
        sleep_time : float, optional
            sleep time, by default 20.0
        refresh_time : float, optional
            refresh time, by default 5.0

        Returns
        -------
        str
            simulation status

        Raises
        ------
        NotImplementedError
            platform not be implemented
        """

        if self.platform == "ubuntu":

            proc = subprocess.Popen(command, shell=True)

            start_time = time.perf_counter()
            time.sleep(sleep_time)
            while True:

                time.sleep(
                    refresh_time - \
                        ((time.perf_counter() - start_time) % refresh_time)
                )
                end_time = time.perf_counter()
                #
                try:
                    file = open(self.job_name + ".msg")
                    word1 = "THE ANALYSIS HAS BEEN COMPLETED"
                    if word1 in file.read():
                        proc.kill()
                        self.kill_abaqus_process()
                        self.flag = "finished"
                        break
                except FileNotFoundError:
                    print(
                        "abaqus license is not enough,"
                        "waiting for license authorization"
                    )
                # over time killing
                if max_time is not inf:
                    if (end_time - start_time) > max_time:
                        proc.kill()
                        self.kill_abaqus_process()
                        print("overtime kill")
                        self.flag = "killed"
                        break
                # kill the process if error in mesh
                try:
                    file = open(self.job_name + ".dat")
                    word1 = "ERROR"
                    if word1 in file.read():
                        proc.kill()
                        self.kill_abaqus_process()
                        self.flag = "error"
                        break
                except FileNotFoundError:
                    print("under modeling and meshing")
            print(f"simulation time :{(end_time - start_time):2f} s")
            # remove files that influence the simulation process
            self.remove_files(directory=os.getcwd())
        elif self.platform == "cluster" or self.platform == "windows":
            # count time for simulation
            start_time = time.perf_counter()
            os.system(command)
            self.flag = "finished"
            end_time = time.perf_counter()
            print(f"simulation time :{(end_time - start_time):2f} s")

        else:
            raise NotImplementedError("platform not be implemented")

        return self.flag

    def make_new_script(
        self,
        file_name: str,
        status: str = "simulation",
    ) -> None:
        """make a small python script for running abaqus

        Parameters
        ----------
        file_name : str
            file name
        status : str, optional
            status of  simulation, by default "simulation"

        """
        # if the platform is windows, change the path to posix
        if self.sim_info["platform"] == "windows":
            path_temp = Path(self.folder_info["script_path"]).as_posix()
            self.folder_info["script_path"] = path_temp  # for flake8

        if status == "simulation":
            with open(file_name, "w") as file:
                file.write("import os \n")
                file.write("import sys \n")
                file.write("import json \n")
                file.write(
                    "sys.path.extend(['"
                    + str(self.folder_info["script_path"])
                    + "']) \n"
                )
                file.write(
                    "from "
                    + str(self.folder_info["sim_path"])
                    + " import "
                    + str(self.folder_info["sim_script"])
                    + "\n"
                )
                line = "file = '" + str(self.abaqus_input_file) + "' \n"
                file.write(line)
                file.write("with open(file, 'r') as f:\n")
                file.write("	dict = json.load(f)\n")
                file.write(str(self.folder_info["sim_script"]) + "(dict)\n")
            file.close()
        elif status == "post_process":
            with open(file_name, "w") as file:
                file.write("import os\n")
                file.write("import sys\n")
                file.write(
                    "sys.path.extend(['"
                    + str(self.folder_info["script_path"])
                    + "']) \n"
                )
                file.write(
                    "from "
                    + str(self.folder_info["post_path"])
                    + " import "
                    + str(self.folder_info["post_script"])
                    + "\n"
                )
                file.write(
                    str(self.folder_info["post_script"])
                    + "('"
                    + str(self.job_name)
                    + "')\n"
                )
            file.close()
        else:
            raise KeyError("provide correct status of simulation \n")

    @staticmethod
    def kill_abaqus_process() -> None:
        """kill abaqus simulation process"""
        standard_solver = "pkill standard"
        os.system(standard_solver)
        ABQcaeK = "pkill ABQcaeK"
        os.system(ABQcaeK)
        SMAPython = "pkill SMAPython"
        os.system(SMAPython)

    @staticmethod
    def remove_files(
        directory: str,
        file_types: list = [
            ".log",
            ".lck",
            ".SMABulk",
            ".rec",
            ".SMAFocus",
            ".exception",
            ".simlog",
            ".023",
            ".exception",
        ],
    ) -> None:
        """remove file

        Parameters
        ----------
        directory : str
            target folder
        file_type : str
            file name
        """
        # get all files in this folder
        all_files = os.listdir(directory)
        for target_file in file_types:
            # get the target file names
            filtered_files = [
                file for file in all_files if file.endswith(target_file)
            ]

            # remove the target files is existed
            for file in filtered_files:
                path_to_file = os.path.join(directory, file)
                if os.path.exists(path_to_file):
                    os.remove(path_to_file)
