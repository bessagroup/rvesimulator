#                                                                       Modules
# =============================================================================
# Standard
import json
import os

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
        dict that constains
    file_name : str
        file name of json file
    """

    with open(file_name, "w") as fp:
        json.dump(sim_info, fp)


def create_working_folder(
    self,
    folder_index: int = None,
    sub_folder_index: int = None,
    third_folder_index: int = None,
) -> None:
    """create working folder

    Parameters
    ----------
    folder_index : int, optional
        folder index, by default None
    sub_folder_index : int, optional
        second folder index , by default None
    third_folder_index : int, optional
        third folder index, by default None

    """
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
