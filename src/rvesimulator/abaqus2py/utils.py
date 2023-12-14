"""
Functions for general use for abaqus2py module
"""
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
    except OSError:
        print(f"Directory {dir_name} can not be created")

    return path


def write_json(sim_info: dict, file_name: str) -> None:
    """write json file for abaqus

    Parameters
    ----------
    sim_info : dict
        dict that contains all simulation information
    file_name : str
        file name of json file
    """

    with open(file_name, "w") as fp:
        json.dump(sim_info, fp)


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
