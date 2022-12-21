import json
import os


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


def make_new_script(
    new_file_name: str,
    folder_info: dict,
    status: str = "simulation",
    sim_info_name: str = "sim_info.json",
    job_name: str = "Job-1",
) -> None:
    """
    Parameters
    ----------
    job_name
    status
    new_file_name
    sim_info_name
    folder_info

    Returns
    -------

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
                str(folder_info["post_script"]) + "('" + str(job_name) + "')\n"
            )
        file.close()


def kill_abaqus_process():
    """_summary_"""
    name_1 = "pkill standard"
    aa = os.system(name_1)
    name_2 = "pkill ABQcaeK"
    bb = os.system(name_2)
    name_3 = "pkill SMAPython"
    cc = os.system(name_3)
    print(aa + bb + cc)


def print_banner(message: str, sign="#", length=50) -> None:
    """_summary_

    Parameters
    ----------
    message : str
        _description_
    sign : str, optional
        _description_, by default "#"
    length : int, optional
        _description_, by default 50
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


def remove_files(directory, file_type: str):
    files_in_directory = os.listdir(directory)
    filtered_files = [
        file for file in files_in_directory if file.endswith(file_type)
    ]
    for file in filtered_files:
        path_to_file = os.path.join(directory, file)
        if os.path.exists(path_to_file):
            os.remove(path_to_file)
