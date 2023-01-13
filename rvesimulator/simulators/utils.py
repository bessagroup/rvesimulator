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


def create_dir(current_folder: str, dirname: str) -> str:
    """
    Make directories and return the path of the directory
    Parameters
    ----------
    current_folder (str): cwd
    dirname (str): folder that needs to be created

    Returns
    -------
    directory path that is created
    """

    path = os.path.join(current_folder, dirname)
    try:
        os.makedirs(path, exist_ok=True)

    except OSError as error:
        print(f"Directory {dirname} can not be created")

    return path


def write_json(sim_info: dict, filename: str) -> None:
    """

    Parameters
    ----------
    sim_info: a dict that contains the information for simulation
    filename: a string of the new .py file

    Returns
    -------
    """
    with open(filename, "w") as fp:
        json.dump(sim_info, fp)
