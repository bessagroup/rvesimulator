"""Abstract class for a FEM simulator and assert inputs
"""
#                                                                       Modules
# =============================================================================

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================

from abc import ABC, abstractmethod


class Simulator(ABC):
    """Base class for a FEM simulator"""

    @abstractmethod
    def pre_process(self,
                    py_script: str,
                    py_func: str) -> None:
        """Function that calls the FEM simulator the pre-processing"""

        raise NotImplementedError("should be implemented in subclass")
    @abstractmethod
    def submit_job(self,
                   num_cpu: int) -> None:
        """Function that calls the FEM simulator to submit the job"""

        raise NotImplementedError("should be implemented in subclass")

    @abstractmethod
    def post_process(self,
                     post_py_script: str,
                     post_py_func: str,
                     delete_odb: bool = True) -> None:
        """Function that handles the post-processing"""

        raise NotImplementedError("should be implemented in subclass")

    def run(self,
            py_script: str,
            py_func: str,
            post_py_script: str = None,
            post_py_func: str = None,
            num_cpu: int = 1,
            delete_odb: bool = True) -> None:
        """run the simulation in one shot
        """
        if post_py_script is None and post_py_func is None:

            self.pre_process(py_script=py_script, py_func=py_func)
            self.post_process(post_py_script=post_py_script,
                              post_py_func=post_py_func,
                              delete_odb=delete_odb)
        else:
            self.pre_process(py_script=py_script, py_func=py_func)
            self.submit_job(num_cpu=num_cpu)
            self.post_process(post_py_script=post_py_script,
                              post_py_func=post_py_func,
                              delete_odb=delete_odb)


class AssertInputs:

    @classmethod
    def is_mwd_in_folder_info(cls, folder_info: dict) -> None:
        """assert main_work_directory in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "main_dir" in folder_info.keys()
        ), "main_dir should in folder_info dict"

    @classmethod
    def is_script_path_in_folder_info(cls, folder_info: dict) -> None:
        """assert script_path in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "script_path" in folder_info.keys()
        ), "script_path should in folder_info dict"



    @classmethod
    def is_sim_script_in_folder_info(cls, folder_info: dict) -> None:
        """assert sim_script in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "sim_script" in folder_info.keys()
        ), "sim_script should in folder_info dict"

    @classmethod
    def is_post_script_in_folder_info(cls, folder_info: dict) -> None:
        """assert post script in the folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "post_script" in folder_info.keys()
        ), "post_script should in folder_info dict"

    @classmethod
    def is_job_name_in_sim_info(cls, sim_info: dict) -> None:
        """assert job name in the folder dict

        Parameters
        ----------
        sim_info : dict
            dict that contains all folder information
        """
        assert (
            "job_name" in sim_info.keys()
        ), "job_name should in folder_info dict"
