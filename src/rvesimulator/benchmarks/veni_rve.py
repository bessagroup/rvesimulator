#                                                                       Modules
# =============================================================================
# Standard
import os
from typing import Any

import rvesimulator
from rvesimulator.additions.microstructure_wrapper import CircleMircoStructure

from .shared_functionalities import SimulationBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================

# =============================================================================
class BothVPLeonovRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.veni_rve",
            "sim_script": "BothVPLeonovRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "/benchmark_abaqus_scripts/vp_leonov_model.f"
        self.update_sim_info()

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            485.77,
            0.45,
            2.47e8,
            4.12e-28,
            7.51e6,
            2.38e-2,
            4.99e1,
            1.00e-2,
            1.03e-1,
            296.0,
            0.1,
        ],
        paras_pe: list = [
            400.14,
            0.45,
            2.34e8,
            4.63e-26,
            6.38e6,
            9.7e-2,
            1.63e1,
            2.56e1,
            1.98e0,
            296.0,
            0.1,
        ],
        mesh_partition: int = 50,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 10.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )

        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform, }

        self.sim_info = {
            "job_name": "veni_vp_leonov_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
            "subroutine_path": self.subroutine_path
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

# =============================================================================
class BothVPLeonovPathLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
            "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.veni_rve",
            "sim_script": "BothVPLeonovPathLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "/benchmark_abaqus_scripts/vp_leonov_model.f"
        self.update_sim_info()

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            485.77,
            0.45,
            2.47e8,
            4.12e-28,
            7.51e6,
            2.38e-2,
            4.99e1,
            1.00e-2,
            1.03e-1,
            296.0,
            0.1,
        ],
        paras_pe: list = [
            400.14,
            0.45,
            2.34e8,
            4.63e-26,
            6.38e6,
            9.7e-2,
            1.63e1,
            2.56e1,
            1.98e0,
            296.0,
            0.1,
        ],
        mesh_partition: int = 50,
        strain: list = [0.05, 0.05, 0.05],
        strain_amplitude: list = None,
        num_steps: int = 1000,
        simulation_time: float = 10.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )

        # simulation paras
        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform, }

        self.sim_info = {
            "job_name": "veni_vp_leonov_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "strain_amplitude": strain_amplitude,
            "num_cpu": num_cpu,
            "platform": platform,
            "subroutine_path": self.subroutine_path
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)


# =============================================================================
class BothVEVPLeonovRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
            "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.veni_rve",
            "sim_script": "BothVEVPLeonovRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "/benchmark_abaqus_scripts/vevp_leonov_model.f"
        self.update_sim_info()

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            748.234,
            0.45,
            296.0,
            0.1,
            1.53e-02,
            6.608e-01,
            7.33e02,
            2.51e08,
            3.30e-26,
            7.42e06,
            2.67e-02,
            4.826e1,
            1.520,
            1.303,
        ],
        paras_pe: list = [
            611.6,
            0.45,
            296.0,
            0.1,
            4.40e-01,
            5.1e-03,
            772.1,
            2.29e08,
            9.94e-26,
            5.01e06,
            3.234e-01,
            15.88,
            13.52,
            0.043,
        ],
        mesh_partition: int = 50,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 1000,
        simulation_time: float = 10.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )
        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform, }

        self.sim_info = {
            "job_name": "veni_ve_vp_leonov_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
            "subroutine_path": self.subroutine_path
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

# =============================================================================
class BothVEVPLeonovPathLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.veni_rve",
            "sim_script": "BothVEVPLeonovPathLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "/benchmark_abaqus_scripts/vevp_leonov_model.f"
        self.update_sim_info()

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            748.234,
            0.45,
            296.0,
            0.1,
            1.53e-02,
            6.608e-01,
            7.33e02,
            2.51e08,
            3.30e-26,
            7.42e06,
            2.67e-02,
            4.826e1,
            1.520,
            1.303,
        ],
        paras_pe: list = [
            611.6,
            0.45,
            296.0,
            0.1,
            4.40e-01,
            5.1e-03,
            772.1,
            2.29e08,
            9.94e-26,
            5.01e06,
            3.234e-01,
            15.88,
            13.52,
            0.043,
        ],
        mesh_partition: int = 50,
        strain: list = [0.05, 0.05, 0.05],
        strain_amplitude: list = None,
        num_steps: int = 1000,
        simulation_time: float = 10.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )

        # simulation paras
        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform, }

        self.sim_info = {
            "job_name": "veni_vevp_leonov_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "strain_amplitude": strain_amplitude,
            "num_cpu": num_cpu,
            "platform": platform,
            "subroutine_path": self.subroutine_path
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)


# =============================================================================
class VEVP_VPLeonovRegularLoads(SimulationBase):
    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                "/scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.veni_rve",
            "sim_script": "VEVP_VPLeonovRegularLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "/benchmark_abaqus_scripts/veni_mix_model.f"
        self.update_sim_info()

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            485.77,
            0.45,
            2.47e8,
            4.12e-28,
            7.51e6,
            2.38e-2,
            4.99e1,
            1.00e-2,
            1.03e-1,
            296.0,
            0.1,
        ],
        paras_pe: list = [
            601.6,
            0.45,
            298.0,
            0.1,
            0.4004,
            0.0006281,
            709.0,
            233000000.0,
            4.0209e-25,
            5495000.0,
            0.3448,
            21.0,
            13.32,
            0.1],
        mesh_partition: int = 50,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 1000,
        simulation_time: float = 10.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )
        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform,
            }

        self.sim_info = {
            "job_name": "veni_mix_leonov_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "num_cpu": num_cpu,
            "platform": platform,
            "subroutine_path": self.subroutine_path
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)


class VEVP_VPLeonovPathLoads(SimulationBase):

    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""
        self.main_folder = os.getcwd()
        self.folder_info = {
            "main_work_directory": os.path.join(os.getcwd(), "Data"),
            "script_path": os.path.dirname(rvesimulator.__file__) + \
                        "scriptbase",
            "current_work_directory": "point_1",
            "sim_path": "benchmark_abaqus_scripts.veni_rve",
            "sim_script": "VEVP_VPLeonovPathLoads",
            "post_path": "basic_analysis_scripts.post_process",
            "post_script": "PostProcess2D",
        }
        self.subroutine_path = self.folder_info["script_path"] + \
            "/benchmark_abaqus_scripts/veni_mix_model.f"
        self.update_sim_info()

    def update_sim_info(
        self,
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        paras_pp: list = [
            485.77,
            0.45,
            2.47e8,
            4.12e-28,
            7.51e6,
            2.38e-2,
            4.99e1,
            1.00e-2,
            1.03e-1,
            296.0,
            0.1],
        paras_pe: list = [
            6.016e+02,
            0.45,
            296.0,
            0.1,
            4.014e-01,
            9.35e-04,
            7.007e+02,
            2.329e+08,
            3.32e-25,
            5.41e+06,
            2.903e-01,
            2.53e+01,
            1.378e+01,
            8.565e-05,
        ],
        mesh_partition: int = 50,
        strain: list = [0.05, 0.05, 0.05],
        strain_amplitude: list = None,
        num_steps: int = 1000,
        simulation_time: float = 10.0,
        num_cpu: int = 1,
        platform: str = "ubuntu",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        # get the micro_structure information
        mirco_structure_generator = CircleMircoStructure()
        microstructure_info, vol_frac = mirco_structure_generator(
            size=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
            seed=seed,
        )

        # simulation paras
        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
            "vol_frac": vol_frac,
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "num_cpu": num_cpu,
            "platform": platform,
            }

        self.sim_info = {
            "job_name": "veni_mix_leonov_rve",
            "location_information": microstructure_info[
                "location_information"
            ],
            "radius_mu": microstructure_info["radius_mu"],
            "radius_std": microstructure_info["radius_std"],
            "len_start": microstructure_info["len_start"],
            "len_end": microstructure_info["len_end"],
            "wid_start": microstructure_info["wid_start"],
            "wid_end": microstructure_info["wid_end"],
            "paras_pp": paras_pp,
            "paras_pe": paras_pe,
            "mesh_partition": mesh_partition,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "strain": strain,
            "strain_amplitude": strain_amplitude,
            "num_cpu": num_cpu,
            "platform": platform,
            "subroutine_path": self.subroutine_path
        }

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)
