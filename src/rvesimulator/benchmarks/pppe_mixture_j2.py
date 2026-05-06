"""
CDDM RVE case
"""
#                                                                       Modules
# =============================================================================
# Standard
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
# local
import rvesimulator
from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
from rvesimulator.microstructure import MicrostructureGenerator, CircleParticles, AmorphousParticles

from .py3rve_base import Py3RVEBase

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "alpha"
# =============================================================================


class J2StrucMesh2DRVE(Py3RVEBase):
    """Interface between python and abaqus of the StrucMesh2D case

    Parameters
    ----------
    SimulationBase : class
        base class for simulation
    """

    def __init__(self) -> None:
        """Interface between python and abaqus of the Hollow plate case"""

        logging.basicConfig(level=logging.INFO, filename="StrucMesh.log")
        self.logger = logging.getLogger("abaqus_simulation")

        self.main_folder = Path.cwd()
        self.folder_info = {
            "main_dir": Path(self.main_folder, str("Data")),
            "script_path": Path(rvesimulator.__file__).parent.as_posix() +
            "/scriptbase",
            "current_dir": "point_1",
            "sim_script": "structural_mesh_scripts.structural_mesh_rve_script",
            "sim_func": "simulation_script",
            "post_script": "basic_analysis_scripts.post_process",
            "post_func": "PostProcess2D",
        }

    def update_sim_info(
        self,
        microstructure_descriptor: MicrostructureGenerator,
        youngs_modulus_matrix: float =  664.03,
        poisson_ratio_matrix: float = 0.45,
        youngs_modulus_fiber: float = 797.97,
        poisson_ratio_fiber: float = 0.42,
        mesh_partition: int = 30,
        strain: List = [0.5, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        hardening_table_fiber: List =[
            [ 4.304324, 0.000000],
            [11.230847, 0.003015],
            [14.370436, 0.006030],
            [16.296567, 0.009045],
            [17.715528, 0.012060],
            [18.874447, 0.015075],
            [19.849797, 0.018090],
            [20.674973, 0.021106],
            [21.366450, 0.024121],
            [21.925699, 0.027136],
            [22.374812, 0.030151],
            [22.715125, 0.033166],
            [22.956229, 0.036181],
            [23.130807, 0.039196],
            [23.248277, 0.042211],
            [23.331979, 0.045226],
            [23.392533, 0.048241],
            [23.437206, 0.051256],
            [23.476556, 0.054271],
            [23.508875, 0.057286],
            [23.539323, 0.060302],
            [23.570658, 0.063317],
            [23.606091, 0.066332],
            [23.645995, 0.069347],
            [23.685668, 0.072362],
            [23.724133, 0.075377],
            [23.766465, 0.078392],
            [23.811855, 0.081407],
            [23.853417, 0.084422],
            [23.899639, 0.087437],
            [23.951407, 0.090452],
            [24.003825, 0.093467],
            [24.056219, 0.096482],
            [24.107561, 0.099497],
            [24.159203, 0.102513],
            [24.205738, 0.105528],
            [24.256461, 0.108543],
            [24.315292, 0.111558],
            [24.373447, 0.114573],
            [24.434232, 0.117588],
            [24.491735, 0.120603],
            [24.545158, 0.123618],
            [24.600396, 0.126633],
            [24.664226, 0.129648],
            [24.732367, 0.132663],
            [24.790409, 0.135678],
            [24.850739, 0.138693],
            [24.913935, 0.141709],
            [24.975599, 0.144724],
            [25.032830, 0.147739],
            [25.097469, 0.150754],
            [25.165154, 0.153769],
            [25.232838, 0.156784],
            [25.300523, 0.159799],
            [25.368207, 0.162814],
            [25.435892, 0.165829],
            [25.503576, 0.168844],
            [25.571261, 0.171859],
            [25.638986, 0.174874],
            [25.706750, 0.177889],
            [25.774419, 0.180905],
            [25.841979, 0.183920],
            [25.909540, 0.186935],
            [25.977100, 0.189950],
            [26.045303, 0.192965],
            [26.113592, 0.195980],
            [26.181881, 0.198995],
            [26.250171, 0.202010],
            [26.318460, 0.205025],
            [26.386753, 0.208040],
            [26.455062, 0.211055],
            [26.523371, 0.214070],
            [26.591680, 0.217085],
            [26.659988, 0.220101],
            [26.728297, 0.223116],
            [26.796606, 0.226131],
            [26.864915, 0.229146],
            [26.933223, 0.232161],
            [27.001532, 0.235176],
            [27.069841, 0.238191],
            [27.138150, 0.241206],
            [27.206459, 0.244221],
            [27.274767, 0.247236],
            [27.343076, 0.250251],
            [27.411385, 0.253266],
            [27.479694, 0.256281],
            [27.548002, 0.259296],
            [27.616311, 0.262312],
            [27.684620, 0.265327],
            [27.752929, 0.268342],
            [27.821237, 0.271357],
            [27.889546, 0.274372],
            [27.957855, 0.277387],
            [28.026164, 0.280402],
            [28.094472, 0.283417],
            [28.162781, 0.286432],
            [28.231090, 0.289447],
            [28.299399, 0.292462],
            [28.367708, 0.295477],
            [28.436016, 0.298492],
            [28.504325, 0.301508],
            [28.572634, 0.304523],
            [28.640943, 0.307538],
            [28.709251, 0.310553],
            [28.777796, 0.313568],
            [28.846460, 0.316583],
            [28.915124, 0.319598],
            [28.983788, 0.322613],
            [29.052452, 0.325628],
            [29.121116, 0.328643],
            [29.189780, 0.331658],
            [29.258444, 0.334673],
            [29.327109, 0.337688],
            [29.395773, 0.340704],
            [29.464437, 0.343719],
            [29.533101, 0.346734],
            [29.601765, 0.349749],
            [29.670429, 0.352764],
            [29.739093, 0.355779],
            [29.807758, 0.358794],
            [29.876422, 0.361809],
            [29.945086, 0.364824],
            [30.013750, 0.367839],
            [30.082414, 0.370854],
            [30.151078, 0.373869],
            [30.219742, 0.376884],
            [30.288406, 0.379900],
            [30.357071, 0.382915],
            [30.425735, 0.385930],
            [30.494399, 0.388945],
            [30.563063, 0.391960],
            [30.631727, 0.394975],
            [30.700391, 0.397990],
            [30.769055, 0.401005],
            [30.837720, 0.404020],
            [30.906384, 0.407035],
            [30.975048, 0.410050],
            [31.043712, 0.413065],
            [31.112376, 0.416080],
            [31.181040, 0.419095],
            [31.249704, 0.422111],
            [31.318368, 0.425126],
            [31.387033, 0.428141],
            [31.455697, 0.431156],
            [31.524361, 0.434171],
            [31.593025, 0.437186],
            [31.661689, 0.440201],
            [31.730353, 0.443216],
            [31.799017, 0.446231],
            [31.867682, 0.449246],
            [31.936346, 0.452261],
            [32.005010, 0.455276],
            [32.073674, 0.458291],
            [32.142338, 0.461307],
            [32.211002, 0.464322],
            [32.279666, 0.467337],
            [32.348330, 0.470352],
            [32.416995, 0.473367],
            [32.485659, 0.476382],
            [32.554323, 0.479397],
            [32.622987, 0.482412],
            [32.691651, 0.485427],
            [32.760315, 0.488442],
            [32.828979, 0.491457],
            [32.897644, 0.494472],
            [32.966308, 0.497487],
            [33.034972, 0.500503],
            [33.103636, 0.503518],
            [33.172300, 0.506533],
            [33.240964, 0.509548],
            [33.309628, 0.512563],
            [33.378292, 0.515578],
            [33.446957, 0.518593],
            [33.515621, 0.521608],
            [33.584285, 0.524623],
            [33.652949, 0.527638],
            [33.721613, 0.530653],
            [33.790277, 0.533668],
            [33.858941, 0.536683],
            [33.927605, 0.539698],
            [33.996270, 0.542714],
            [34.064934, 0.545729],
            [34.133598, 0.548744],
            [34.202262, 0.551759],
            [34.270926, 0.554774],
            [34.339590, 0.557789],
            [34.408254, 0.560804],
            [34.476919, 0.563819],
            [34.545583, 0.566834],
            [34.614247, 0.569849],
            [34.682911, 0.572864],
            [34.751575, 0.575879],
            [34.820239, 0.578894],
            [34.888903, 0.581910],
            [34.957568, 0.584925],
            [35.026232, 0.587940],
            [35.094896, 0.590955],
            [35.163560, 0.593970],
            [35.232224, 0.596985],
            [35.300888, 0.60000 ],
        ],
        hardening_table_matrix: list =[
            [3.023663, 0.000000 ],
            [16.379144, 0.020101 ],
            [19.944192, 0.040201 ],
            [21.955358, 0.060302 ],
            [23.223464, 0.080402 ],
            [24.065270, 0.100503 ],
            [24.639390, 0.120603 ],
            [25.044329, 0.140704 ],
            [25.336121, 0.160804 ],
            [25.568512, 0.180905 ],
            [25.743041, 0.201005 ],
            [25.870881, 0.221106 ],
            [25.983502, 0.241206 ],
            [26.049062, 0.261307 ],
            [26.094383, 0.281407 ],
            [26.090894, 0.301508 ],
            [26.007020, 0.321608 ],
            [25.868177, 0.341709 ],
            [25.694250, 0.361809 ],
            [25.500848, 0.381910 ],
            [25.302199, 0.402010 ],
            [25.151576, 0.422111 ],
            [25.060858, 0.442211 ],
            [25.014115, 0.462312 ],
            [25.045785, 0.482412 ],
            [25.162078, 0.502513 ],
            [25.331807, 0.522613 ],
            [25.541338, 0.542714 ],
            [25.795063, 0.562814 ],
            [26.089967, 0.582915 ],
            [26.444631, 0.603015 ],
            [26.797472, 0.623116 ],
            [27.224121, 0.643216 ],
            [27.695375, 0.663317 ],
            [28.197757, 0.683417 ],
            [28.719438, 0.703518 ],
            [29.300995, 0.723618 ],
            [29.901527, 0.743719 ],
            [30.529320, 0.763819 ],
            [31.144322, 0.783920 ],
            [31.785368, 0.804020 ],
            [32.427874, 0.824121 ],
            [33.077093, 0.844221 ],
            [33.773771, 0.864322 ],
            [34.477562, 0.884422 ],
            [35.176534, 0.904523 ],
            [35.891482, 0.924623 ],
            [36.662110, 0.944724 ],
            [37.454270, 0.964824 ],
            [38.198240, 0.984925 ],
            [39.114768, 1.005025 ],
            [40.467594, 1.025126 ],
            [41.820420, 1.045226 ],
            [43.226816, 1.065327 ],
            [44.660296, 1.085427 ],
            [46.094507, 1.105528 ],
            [47.529446, 1.125628 ],
            [48.965118, 1.145729 ],
            [50.401520, 1.165829 ],
            [51.838651, 1.185930 ],
            [53.276512, 1.206030 ],
            [54.715102, 1.226131 ],
            [56.154423, 1.246231 ],
            [57.594475, 1.266332 ],
            [59.054640, 1.286432 ],
            [60.537391, 1.306533 ],
            [62.040743, 1.326633 ],
            [63.564763, 1.346734 ],
            [65.109453, 1.366834 ],
            [66.674745, 1.386935 ],
            [68.260727, 1.407035 ],
            [69.867353, 1.427136 ],
            [71.494589, 1.447236 ],
            [73.142526, 1.467337 ],
            [74.811092, 1.487437 ],
            [76.500276, 1.507538 ],
            [78.210164, 1.527638 ],
            [79.940666, 1.547739 ],
            [81.691801, 1.567839 ],
            [83.463640, 1.587940 ],
            [85.256081, 1.608040 ],
            [87.069177, 1.628141 ],
            [88.902952, 1.648241 ],
            [90.757331, 1.668342 ],
            [92.632390, 1.688442 ],
            [94.528105, 1.708543 ],
            [96.444422, 1.728643 ],
            [98.381440, 1.748744 ],
            [100.339092, 1.768844],
            [102.317358, 1.788945],
            [104.316329, 1.809045],
            [106.342713, 1.829146],
            [108.583239, 1.849246],
            [110.854437, 1.869347],
            [113.156211, 1.889447],
            [115.488612, 1.909548],
            [117.851666, 1.929648],
            [120.245285, 1.949749],
            [122.669537, 1.969849],
            [125.124444, 1.989950],
            [127.609917, 2.010050],
            [130.126040, 2.030151],
            [132.672775, 2.050251],
            [135.250102, 2.070352],
            [138.067949, 2.090452],
            [141.094091, 2.110553],
            [144.165939, 2.130653],
            [147.371566, 2.150754],
            [150.869077, 2.170854],
            [154.425093, 2.190955],
            [158.039727, 2.211055],
            [161.712853, 2.231156],
            [165.444553, 2.251256],
            [169.234808, 2.271357],
            [173.083586, 2.291457],
            [176.990941, 2.311558],
            [180.956856, 2.331658],
            [184.981256, 2.351759],
            [189.064244, 2.371859],
            [193.205808, 2.391960],
            [197.405859, 2.412060],
            [201.664533, 2.432161],
            [205.981686, 2.452261],
            [210.357404, 2.472362],
            [214.791704, 2.492462],
            [219.284537, 2.512563],
            [223.835884, 2.532663],
            [228.445851, 2.552764],
            [233.114307, 2.572864],
            [237.841315, 2.592965],
            [242.626896, 2.613065],
            [247.470993, 2.633166],
            [252.373693, 2.653266],
            [257.334901, 2.673367],
            [262.354634, 2.693467],
            [267.432960, 2.713568],
            [272.569842, 2.733668],
            [277.765215, 2.753769],
            [283.019206, 2.773869],
            [288.331710, 2.793970],
            [293.702738, 2.814070],
            [299.132369, 2.834171],
            [304.620495, 2.854271],
            [310.167205, 2.874372],
            [315.772469, 2.894472],
            [321.436258, 2.914573],
            [327.158582, 2.934673],
            [332.939499, 2.954774],
            [338.778919, 2.974874],
            [344.676939, 2.994975],
            [350.633478, 3.015075],
            [356.648525, 3.035176],
            [362.722202, 3.055276],
            [368.854372, 3.075377],
            [375.045083, 3.095477],
            [381.294389, 3.115578],
            [387.602234, 3.135678],
            [393.968579, 3.155779],
            [400.393534, 3.175879],
            [406.876985, 3.195980],
            [413.419023, 3.216080],
            [420.019619, 3.236181],
            [426.678703, 3.256281],
            [433.396386, 3.276382],
            [440.172608, 3.296482],
            [447.007355, 3.316583],
            [453.900679, 3.336683],
            [460.852562, 3.356784],
            [467.862931, 3.376884],
            [474.931938, 3.396985],
            [482.059423, 3.417085],
            [489.245460, 3.437186],
            [496.490111, 3.457286],
            [503.793253, 3.477387],
            [511.154943, 3.497487],
            [518.575200, 3.517588],
            [526.053992, 3.537688],
            [533.591336, 3.557789],
            [541.187264, 3.577889],
            [548.841664, 3.597990],
            [556.554696, 3.618090],
            [564.326236, 3.638191],
            [572.156284, 3.658291],
            [580.044957, 3.678392],
            [587.992155, 3.698492],
            [595.997865, 3.718593],
            [604.062170, 3.738693],
            [612.184991, 3.758794],
            [620.366346, 3.778894],
            [628.606331, 3.798995],
            [636.904780, 3.819095],
            [645.261805, 3.839196],
            [653.677401, 3.859296],
            [662.151499, 3.879397],
            [670.684182, 3.899497],
            [679.275409, 3.919598],
            [687.925163, 3.939698],
            [696.633499, 3.959799],
            [705.400372, 3.979900],
            [714.225708, 4.00000 ],
        ],
        element_type: str = "CPE4R",
        solver_type: str = "nonlinear",
        seed: Any = None,
        print_info: bool = False,
    ) -> None:
        """regular rve for Ti6Al4V

        Parameters
        ----------
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
        hardening_law_fiber : Any, optional
            hardening law for the fiber material,
            by default LinearHardeningLaw()
        haderning_law_matrix : Any, optional
            hardening law for the matrix material,
            by default LinearHardeningLaw()
        seed : Any, optional
            seed number, by default None
        print_info : bool, optional
            print simulation information or not, by default False
        """

        # get simulation information
        self.microstructure_descriptor = microstructure_descriptor
        self.youngs_modulus_matrix = youngs_modulus_matrix
        self.poisson_ratio_matrix = poisson_ratio_matrix
        self.youngs_modulus_fiber = youngs_modulus_fiber
        self.poisson_ratio_fiber = poisson_ratio_fiber
        self.mesh_partition = mesh_partition
        self.strain = strain
        self.num_steps = num_steps
        self.simulation_time = simulation_time
        self.num_cpu = num_cpu
        self.seed = seed
        self.element_type = element_type
        self.solver_type = solver_type
        # create hardening table for matrix
        self.hardening_table_matrix = np.array(
            hardening_table_matrix).reshape(-1, 2).T.tolist()
        # create hardening table for fiber
        self.hardening_table_fiber = np.array(
            hardening_table_fiber).reshape(-1, 2).T.tolist()

        self.sim_paras = {
            "youngs_modulus_matrix": youngs_modulus_matrix,
            "poisson_ratio_matrix": poisson_ratio_matrix,
            "youngs_modulus_fiber": youngs_modulus_fiber,
            "poisson_ratio_fiber": poisson_ratio_fiber,
            "hardening_table_fiber": self.hardening_table_fiber,
            "hardening_table_matrix": self.hardening_table_matrix,
            "mesh_partition": mesh_partition,
            "strain": strain,
            "num_steps": num_steps,
            "simulation_time": simulation_time,
            "solver_type": solver_type,
            "element_type": element_type,
            "num_cpu": num_cpu}

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""

        self.sim_info = {
            "job_name": "J2StrucMesh",
            "length": self.microstructure_descriptor.length,
            "width": self.microstructure_descriptor.width,
            "youngs_modulus_matrix": self.youngs_modulus_matrix,
            "poisson_ratio_matrix": self.poisson_ratio_matrix,
            "youngs_modulus_fiber": self.youngs_modulus_fiber,
            "poisson_ratio_fiber": self.poisson_ratio_fiber,
            "mesh_partition": self.mesh_partition,
            "hardening_table_fiber": self.hardening_table_fiber,
            "hardening_table_matrix": self.hardening_table_matrix,
            "num_steps": self.num_steps,
            "simulation_time": self.simulation_time,
            "strain": self.strain,
            "num_cpu": self.num_cpu,
            "element_type": self.element_type,
            "solver_type": self.solver_type,
        }

    def run_simulation(
        self,
        sample: Dict = None,
        folder_index: int = None,
        delete_odb: bool = True,
    ) -> Dict:
        """run single simulation

        Parameters
        ----------
        sample : dict, optional
            a dict contains the information of design variables
        folder_index : int, optional
            first folder index, by default None
        sub_folder_index : int, optional
            second folder index, by default None
        third_folder_index : int, optional
            third folder index, by default None

        Returns
        -------
        dict
            all the simulation results from abaqus
        """
        # number of samples
        self._create_working_folder(
            folder_index,
        )
        self.logger.info("working folder: {}".format(self.working_folder))
        self.microstructure_descriptor.generate_microstructure(seed=self.seed)
        self.microstructure_descriptor.to_abaqus_format()
        self.microstructure_descriptor.plot_microstructure(save_figure=True,
                                                fig_name="rve_{}.png".
                                                format(self.seed))
        if isinstance(self.microstructure_descriptor, CircleParticles):
            self.microstructure_descriptor.crate_rgmsh(
                num_discrete=self.mesh_partition)
        elif isinstance(self.microstructure_descriptor, AmorphousParticles):
            self.microstructure_descriptor.crate_rgmsh()
        else:
            raise ValueError("microstructure_descriptor is not supported")
        self.microstructure_descriptor.to_crate_format()
        # check the "microstructure.rgmsh.npy" file is in the folder or not
        if not os.path.exists("microstructure.rgmsh.npy"):
            self.logger.info("microstructure.rgmsh.npy is not in the folder")
            print("microstructure.rgmsh.npy is not in the folder")
        else:
            self.logger.info("microstructure.rgmsh.npy is in the folder")
            print("microstructure.rgmsh.npy is in the folder")

        self.vol_frac = self.microstructure_descriptor.vol_frac
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
                          post_py_script=self.folder_info["post_script"],
                          num_cpu=self.num_cpu,
                          delete_odb=delete_odb)
            results = simulator.read_back_results()
            self.logger.info("abaqus simulation finished")
        except FileNotFoundError:
            self.logger.info("abaqus simulation failed")
            results = None
        # get the simulation results back
        end_time = time.time()
        self.logger.info("time used: {} s".format(end_time - start_time))
        self.logger.info("============== End abaqus simulation ============")

        # back to main folder
        os.chdir(self.main_folder)

        return results
