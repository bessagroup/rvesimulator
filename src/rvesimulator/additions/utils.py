"""
Functions for general use for abaqus2py module
"""
#                                                                       Modules
# =============================================================================
# Standard
import json
import os
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


def create_dir(current_folder: Path,
               dir_name: str) -> Path:
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

    path = current_folder / dir_name

    # create the folder if it does not exist
    path.mkdir(parents=True, exist_ok=True)

    return path


def write_json(sim_info: Dict,
               file_name: str) -> None:
    """write json file for abaqus
    """

    with open(file_name, "w") as fp:
        json.dump(sim_info, fp)


def print_banner(message: str,
                 sign="#", length=50) -> None:
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
    # print(sign * length)


def make_new_script(
    file_name: str,
    script_path: str,
    script_name: str,
    func_name: str,
    input_file_name: str,
) -> None:
    """make a small python script for running abaqus

    """

    with open(file_name, "w") as file:
        file.write("import os \n")
        file.write("import sys \n")
        file.write("import json \n")
        file.write(
            "sys.path.extend(['"
            + str(script_path)
            + "']) \n"
        )
        file.write(
            "from "
            + str(script_name)
            + " import "
            + str(func_name)
            + "\n"
        )
        line = "file = '" + str(input_file_name) + "' \n"
        file.write(line)
        file.write("with open(file, 'r') as f:\n")
        file.write("	dict = json.load(f)\n")
        file.write(str(func_name) + "(dict)\n")
    file.close()


def remove_files(
    directory: str,
    file_types: List = [
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


def rve_microstructure_plot(
    fibers: np.ndarray,
    size: float,
    save_fig: bool = False,
    fig_name: str = "rve.png",
    **kwargs,
) -> None:
    """plot the microstructure if needed

    Parameters
    ----------
    fibers : np.ndarray
        fiber locations
    size : float
        size of rve
    save_fig : bool, optional
        save figure, by default False
    fig_name : str, optional
        figure name, by default "rve.png"
    """

    fig, axes = plt.subplots(**kwargs)
    for ii in range(fibers.shape[0]):
        cc = plt.Circle(
            (fibers[ii, 0], fibers[ii, 1]),
            fibers[ii, 2],
            color="#77AADD",
        )
        axes.add_artist(cc)
    axes.set_aspect(1)
    plt.xlim((0.0, size))
    plt.ylim((0.0, size))
    axes.set_yticks([])
    axes.set_xticks([])
    if save_fig is True:
        plt.savefig(fig_name, dpi=300, bbox_inches="tight")
        plt.close()
    plt.show()
