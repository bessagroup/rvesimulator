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


class Simulator:
    """Base class for a FEM simulator"""

    def execute(self) -> None:
        """Function that calls the FEM simulator the pre-processing"""

        raise NotImplementedError("should be implemented in subclass")

    def post_process(self) -> None:
        """Function that handles the post-processing"""

        raise NotImplementedError("should be implemented in subclass")

    def run(self) -> None:
        """run the simulation in one shot
        """
        self.execute()
        self.post_process()


class AssertInputs:
    @classmethod
    def is_inputs_proper_defined(
        cls, folder_info: dict, sim_info: dict
    ) -> None:
        """ assert inputs are properly defined
        """
        cls.is_mwd_in_folder_info(folder_info=folder_info)
        cls.is_script_path_in_folder_info(folder_info=folder_info)
        cls.is_cwd_in_folder_info(folder_info=folder_info)
        cls.is_sim_path_in_folder_info(folder_info=folder_info)
        cls.is_sim_script_in_folder_info(folder_info=folder_info)
        cls.is_post_path_in_folder_info(folder_info=folder_info)
        cls.is_post_script_in_folder_info(folder_info=folder_info)
        cls.is_job_name_in_sim_info(sim_info=sim_info)
        cls.is_platform_in_sim_info(sim_info=sim_info)

    @classmethod
    def is_mwd_in_folder_info(cls, folder_info: dict) -> None:
        """assert main_work_directory in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "main_work_directory" in folder_info.keys()
        ), "main_work_directory should in folder_info dict"

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
    def is_cwd_in_folder_info(cls, folder_info: dict) -> None:
        """assert current_work_directory in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "current_work_directory" in folder_info.keys()
        ), "current_work_directory should in folder_info dict"

    @classmethod
    def is_sim_path_in_folder_info(cls, folder_info: dict) -> None:
        """assert sim_path in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "sim_path" in folder_info.keys()
        ), "sim_path should in folder_info dict"

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
    def is_post_path_in_folder_info(cls, folder_info: dict) -> None:
        """assert post_path in folder_info dict

        Parameters
        ----------
        folder_info : dict
            dict that contains all folder information
        """
        assert (
            "post_path" in folder_info.keys()
        ), "post_path should in folder_info dict"

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

    @classmethod
    def is_platform_in_sim_info(cls, sim_info: dict) -> None:
        """assert platform in the folder dict

        Parameters
        ----------
        sim_info : dict
            dict that contains all folder information
        """
        assert (
            "platform" in sim_info.keys()
        ), "platform should in folder_info dict"
