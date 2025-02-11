"""
Class for replicating the material property linkages from the paper " Bayesian
neural network for uncertainty quantification in data-driven materials modeling"
by Oliver et al. (2021).
"""
#                                                                       Modules
# =============================================================================
# Standard
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, Tuple
import numpy as np
# local
import rvesimulator
from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
from rvesimulator.additions.hardening_law import SwiftHardeningLaw
from rvesimulator.microstructure.circle_particles import CircleParticles
import cratepy
from .py3rve_base import Py3RVEBase
from rvesimulator.rve2crateio.rvesimulator2crateio import rvesimulator2crateIO

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


class StrucPropDNS(Py3RVEBase):
    """Interface between python and Abaqus of material structure property
    linkages


    Parameters
    ----------
    SimulationBase : class
        base class for simulation
    """

    def __init__(self) -> None:
        """Interface between python and abaqus of the  case"""

        logging.basicConfig(level=logging.INFO, filename="StrucProp.log")
        self.logger = logging.getLogger("abaqus_simulation")

        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "structural_mesh_scripts.structure_property_linkage_sve_script",
            "sim_func": "simulation_script",
            "post_script": "structural_mesh_scripts.structure_property_linkage_sve_script",
            "post_func": "post_process",
        }

    def update_sim_info(
        self,
        size: float = 1.0,
        radius_mu: float = 0.03125,
        radius_std: float = 0.003125,
        vol_req: float = 0.30,
        youngs_modulus_matrix: float = 100000.0,
        poisson_ratio_matrix: float = 0.3,
        youngs_modulus_fiber: float = 400000.0,
        poisson_ratio_fiber: float = 0.25,
        mesh_partition: int = 100,
        strain: list = [0.05, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        yield_stress: float = 400,
        param_a: float = 400,
        param_b: float = 0.3,
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        """path dependent rve for cddm

        Parameters
        ----------
        size : float, optional
            size of the rve, by default 0.048
        radius_mu : float, optional
            radius mean of the rve, by default 0.003
        radius_std : float, optional
            radius deviation of the rve, by default 0.0
        vol_req : float, optional
            volume fraction requirement, by default 0.30
        youngs_modulus_matrix : float, optional
            youngs modulus of the matrix material, by default 100.0
        poisson_ratio_matrix : float, optional
            poisson ratio of the matrix material, by default 0.3
        youngs_modulus_fiber : float, optional
            youngs modulus of the fiber material, by default 1.0
        poisson_ratio_fiber : float, optional
            poisson ratio of the fiber material, by default 0.19
        mesh_partition : int, optional
            mesh partition for the edges, by default 30
        strain : list, optional
            applied maximum strain, by default [0.1, 0.0, 0.0]
        num_steps : int, optional
            number of simulation steps, by default 100
        simulation_time : float, optional
            total simulation time, by default 1.0
        num_cpu : int, optional
            number of cpu used for simulation, by default 1
        platform : str, optional
            platform for simulation, by default "ubuntu"
        hardening_law : class, optional
            hardening law for the simulation, by default LinearHardeningLaw()
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get simulation information
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        self.youngs_modulus_matrix = youngs_modulus_matrix
        self.poisson_ratio_matrix = poisson_ratio_matrix
        self.youngs_modulus_fiber = youngs_modulus_fiber
        self.poisson_ratio_fiber = poisson_ratio_fiber
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.hardening_law = SwiftHardeningLaw(yield_stress=yield_stress,
                                               a=param_a,
                                               b=param_b)
        self.seed = seed
        # get hardening law
        self.hardening_table = self.hardening_law.calculate_hardening_table()

        self.sim_params = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "youngs_modulus_matrix": youngs_modulus_matrix,
            "poisson_ratio_matrix": poisson_ratio_matrix,
            "youngs_modulus_fiber": youngs_modulus_fiber,
            "poisson_ratio_fiber": poisson_ratio_fiber,
            "hardening_table_matrix": self.hardening_table,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu, }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_params)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        self.sim_info = {
            "job_name": "StrucProp",
            "location_information": self.microstructure.microstructure_info[
                "location_information"
            ],
            "radius_mu": self.microstructure.microstructure_info["radius_mu"],
            "radius_std": self.microstructure.microstructure_info[
                "radius_std"],
            "len_start": self.microstructure.microstructure_info["len_start"],
            "len_end": self.microstructure.microstructure_info["len_end"],
            "wid_start": self.microstructure.microstructure_info["wid_start"],
            "wid_end": self.microstructure.microstructure_info["wid_end"],
            "youngs_modulus_matrix": self.youngs_modulus_matrix,
            "poisson_ratio_matrix": self.poisson_ratio_matrix,
            "youngs_modulus_fiber": self.youngs_modulus_fiber,
            "poisson_ratio_fiber": self.poisson_ratio_fiber,
            "mesh_partition": self.mesh_partition,
            "hardening_table_matrix": self.hardening_table,
            "num_steps": self.num_steps,
            "simulation_time": self.simulation_time,
            "strain": self.strain,
            "num_cpu": self.num_cpu,
        }

    def run_simulation(
        self,
        sample: Dict = None,
        folder_index: int = 0,
        delete_odb: bool = False,
    ) -> Dict:
        """run single simulation

        Parameters
        ----------
        sample : dict, optional
            a dict contains the information of design variables
        folder_index : int, optional
            first folder index, by default None
        delete_odb : bool, optional
            delete odb file or not, by default False

        Returns
        -------
        Dict
            all the simulation results from abaqus
        """
        # number of samples
        self._create_working_folder(folder_index)

        self.logger.info("working folder: {}".format(self.working_folder))
        # create microstructure
        self.microstructure = CircleParticles(
            length=self.size,
            width=self.size,
            radius_mu=self.radius_mu,
            radius_std=self.radius_std,
            vol_req=self.vol_req,
        )
        self.microstructure.generate_microstructure(seed=self.seed)
        self.microstructure.to_abaqus_format()
        self.microstructure.plot_microstructure(save_figure=True,
                                                fig_name="rve_{}.png".
                                                format(self.seed))
        # generate the structure mesh
        # generate the discrete microstructure
        self.microstructure.crate_rgmsh(num_discrete=self.mesh_partition)
        self.microstructure.to_crate_format()
        # check the "microstructure.rgmsh.npy" file is in the folder or not
        if not os.path.exists("microstructure.rgmsh.npy"):
            self.logger.info("microstructure.rgmsh.npy is not in the folder")
            print("microstructure.rgmsh.npy is not in the folder")
        else:
            self.logger.info("microstructure.rgmsh.npy is in the folder")
            print("microstructure.rgmsh.npy is in the folder")
        self.vol_frac = self.microstructure.vol_frac
        self.logger.info("volume fraction: {}".format(self.vol_frac))
        # update simulation information
        self._get_sim_info()
        # update the geometry info for microstructure
        self._update_sample_info(sample=sample)
        # update logger on samples
        self.logger.info("==============        update info      ============")
        self.logger.info("sample: {}".format(sample))
        # change folder to main folder
        # save microstructure
        # update simulation information
        self.logger.info("============== Start abaqus simulation ============")
        start_time = time.time()
        simulator = AbaqusSimulator(
            sim_info=self.sim_info, folder_info=self.folder_info
        )
        # run abaqus simulation
        try:
            simulator.run(py_func=self.folder_info["sim_func"],
                          py_script=self.folder_info["sim_script"],
                          post_py_func=self.folder_info["post_func"],
                          num_cpu=self.num_cpu,
                          post_py_script=self.folder_info["post_script"],
                          delete_odb=delete_odb)
            # get the simulation results back
            results = simulator.read_back_results()
            self.logger.info("simulation finished")
        except FileNotFoundError:
            self.logger.error("simulation failed")
            results = None
        end_time = time.time()
        self.logger.info("time used: {} s".format(end_time - start_time))
        self.logger.info("============== End abaqus simulation ============")

        # back to main folder
        os.chdir(self.main_folder)

        return results


class StrucPropSVESCA:
    """ Structure property linkage simulation with tension loading using SCA
    """

    def __init__(self,
                 size: float = 1.0,
                 radius_mu: float = 0.03125,
                 radius_std: float = 0.003125,
                 vol_req: float = 0.30,
                 youngs_modulus_matrix: float = 100000.0,
                 poisson_ratio_matrix: float = 0.3,
                 youngs_modulus_fiber: float = 400000.0,
                 poisson_ratio_fiber: float = 0.25,
                 mesh_partition: int = 200,
                 strain: list = [0.05, 0.0, 0.0],
                 num_steps: int = 100,
                 yield_stress: float = 400,
                 param_a: float = 400,
                 param_b: float = 0.3,
                 min_dist_factor: float = 1.3,
                 job_number: int = 1,
                 delete_vtk: bool = True,
                 print_info: bool = False,) -> None:
        """Initialize the simulation

        Parameters
        ----------
        size : float
            size of the microstructure
        radius_mu : float
            mean radius of the particles
        radius_std : float
            standard deviation of the radius of the particles
        vol_req : float
            volume fraction of the particles

        """

        # information of microstructure
        self.size = size
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        self.min_dist_factor = min_dist_factor
        self.mesh_partition = mesh_partition

        # information of material properties
        self.youngs_modulus_matrix = youngs_modulus_matrix
        self.poisson_ratio_matrix = poisson_ratio_matrix
        self.youngs_modulus_fiber = youngs_modulus_fiber
        self.poisson_ratio_fiber = poisson_ratio_fiber
        self.hardening_law = SwiftHardeningLaw(
            a=param_a,
            b=param_b,
            yield_stress=yield_stress)

        # path-dependent loading
        self.strain = strain
        self.num_steps = num_steps

        # print information or not
        self.print_info = print_info
        self.delete_vtk = delete_vtk

        # define the job number
        self.job_number = job_number

    def configure_input_file(self,
                             clusters: list = [2, 1],
                             rve_seed: int = 1,
                             template_path: str = "mat_struc_linkage_input_template.dat",
                             out_file_path: str = "Data") -> Tuple[str, str]:
        """configure the input file for the simulation

        Parameters
        ----------
        template_path : str
            relative path of the template file in the repo
        out_file_path : str
            relative path of the output file in the repo
        """
        # the path of the repo
        file_path = Path(__file__).parents[1].as_posix()
        # get the path of the template file
        absolute_template_path = file_path + "/rve2crateio/" + template_path
        # get the path of current file
        self.main_folder = Path.cwd()
        # create a folder according to out_file_path
        out_folder = Path(self.main_folder, str(out_file_path))
        out_folder.mkdir(parents=True, exist_ok=True)
        # change the directory to the output folder
        os.chdir(out_folder)
        self.cwd = out_folder.as_posix()
        # configure crate adapter
        crate_adapter = rvesimulator2crateIO()
        # read the template file
        try:
            crate_adapter.read_template(absolute_template_path)
        except FileNotFoundError:
            print("Error: template file not found")
            return
        # generate the microstructure
        microstructure = CircleParticles(
            length=self.size,
            width=self.size,
            radius_mu=self.radius_mu,
            radius_std=self.radius_std,
            vol_req=self.vol_req,
            dist_min_factor=self.min_dist_factor,
        )
        # generate microstructure and write to file
        microstructure.generate_microstructure(seed=rve_seed)
        # plot the microstructure in the output folder
        microstructure.plot_microstructure(
            save_figure=True,
            fig_name=f"rve_{rve_seed}.png")
        microstructure.crate_rgmsh(num_discrete=self.mesh_partition)
        microstructure.to_crate_format(
            file_name=f"microstructure_{rve_seed}.rgmsh.npy")
        # configure the microstructure info of the input file
        crate_adapter.change_size_of_rve(self.size)
        crate_adapter.replace_discretization_file(
            file_path=f"microstructure_{rve_seed}.rgmsh.npy")
        # configure the material properties of the input file
        self.hardening_law.calculate_hardening_table()
        table = self.hardening_law.hardening_law_table
        # new table
        new_table = np.zeros((table.shape[1], 2))
        new_table[:, 0] = table[1, :]
        new_table[:, 1] = table[0, :]
        # change the material properties of the matrix
        crate_adapter.change_matrix_material_properties(
            young_modulus=self.youngs_modulus_matrix,
            poisson_ratio=self.poisson_ratio_matrix,
            hardening_law_table=new_table)
        # change the material properties of the fiber
        crate_adapter.change_particle_material_properties(
            young_modulus=self.youngs_modulus_fiber,
            poisson_ratio=self.poisson_ratio_fiber)

        # change the clustering scheme
        crate_adapter.change_clustering_scheme(
            n_matrix_clusters=clusters[0],
            n_particle_clusters=clusters[1])
        # write the input file
        output_file = f"job_{self.job_number}_rve_{rve_seed}_cluster_{clusters[0]}_{clusters[1]}.dat"
        sim_file_path = Path(output_file).absolute()
        crate_adapter.write_file(output_file)
        # print information
        if self.print_info:
            if Path(output_file).exists():
                # check the input file is successfully written or not
                print(f"Input file is written to {output_file}")
                print(f"Path: {Path(output_file).absolute()}")
            else:
                print(f"Error: input file is not written to {output_file}")

        # change the directory back to the main folder
        os.chdir(self.main_folder)

        self.sim_file_path = sim_file_path
        return sim_file_path.as_posix(), self.cwd

    def run_simulation(self
                       ) -> None:
        """run the simulation"""
        cratepy.crate_simulation(arg_input_file_path=self.sim_file_path,
                                 arg_discret_file_dir=self.cwd,
                                 is_null_stdout=False)
