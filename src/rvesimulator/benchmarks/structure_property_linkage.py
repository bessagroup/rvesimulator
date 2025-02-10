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
from typing import Any, Dict

# local
import rvesimulator
from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
from rvesimulator.additions.hardening_law import SwiftHardeningLaw
from rvesimulator.microstructure.circle_particles import CircleParticles

from .py3rve_base import Py3RVEBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================


class StrucPropSVE(Py3RVEBase):
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
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
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
