#                                                                       Modules
# =============================================================================
# Standard
import json
import os

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


def create_dir(current_folder: str, dir_name: str) -> str:
    """create new directory

    Parameters
    ----------
    current_folder : str
        current working folder
    dirname : str
        new folder name

    Returns
    -------
    str
        path of created folder
    """

    path = os.path.join(current_folder, dir_name)
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as error:
        print(f"Directory {dir_name} can not be created")

    return path


def write_json(sim_info: dict, file_name: str) -> None:
    """write json file for abaqus

    Parameters
    ----------
    sim_info : dict
        dict that constains
    file_name : str
        file name of json file
    """

    with open(file_name, "w") as fp:
        json.dump(sim_info, fp)


def create_working_folder(
    self,
    folder_index,
    sub_folder_index,
    third_folder_index,
) -> None:
    if folder_index is None:
        if sub_folder_index is None:
            self.folder_info["current_work_directory"] = "rate_" + str(
                third_folder_index
            )
        else:
            if third_folder_index is None:
                self.folder_info["current_work_directory"] = "point_" + str(
                    sub_folder_index
                )
            else:
                self.folder_info["current_work_directory"] = (
                    "point_"
                    + str(sub_folder_index)
                    + "/case_"
                    + str(third_folder_index)
                )
    else:
        if sub_folder_index is None:
            raise ValueError("provide sub_folder_index")
        elif third_folder_index is None:
            raise ValueError("provide third_folder_index")
        else:
            self.folder_info["current_work_directory"] = (
                "gen_"
                + str(folder_index)
                + "/point_"
                + str(sub_folder_index)
                + "/case_"
                + str(third_folder_index)
            )
    new_path = create_dir(
        current_folder=self.folder_info["main_work_directory"],
        dirname=self.folder_info["current_work_directory"],
    )
    self.working_folder = new_path
    os.chdir(new_path)
