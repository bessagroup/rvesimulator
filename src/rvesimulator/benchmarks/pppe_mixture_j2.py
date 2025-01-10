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
from typing import Any, Dict

import numpy as np
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
        size: float = 0.048,
        radius_mu: float = 0.003,
        radius_std: float = 0.0,
        vol_req: float = 0.30,
        youngs_modulus_matrix: float = 510.36,
        poisson_ratio_matrix: float = 0.45,
        youngs_modulus_fiber: float = 713.13,
        poisson_ratio_fiber: float = 0.45,
        mesh_partition: int = 30,
        strain: list = [0.1, 0.0, 0.0],
        num_steps: int = 100,
        simulation_time: float = 1.0,
        num_cpu: int = 1,
        hardening_law_fiber: Any = SwiftHardeningLaw(),
        hardening_law_matrix: Any = SwiftHardeningLaw(),
        seed: Any = None,
        print_info: bool = False,
        min_dist_factor: float = 1.2,
    ) -> None:
        """regular rve for Ti6Al4V

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
        self.hardening_law_fiber = hardening_law_fiber
        self.hardening_law_matrix = hardening_law_matrix
        self.seed = seed

        # for now hard code the hardening table for fiber and matrix
        pp_table = [6.201145019, 0.0,
                    7.406784123, 0.000222709,
                    7.785127381, 0.000466271,
                    8.138957899, 0.000756893,
                    8.512860278, 0.001007219,
                    8.923006579, 0.001185563,
                    9.30323792, 0.001421559,
                    9.658119036, 0.001706266,
                    10.03280602, 0.001951205,
                    10.41534024, 0.00217981,
                    10.76908969, 0.002463862,
                    11.11595909, 0.002760439,
                    11.42160461, 0.003136838,
                    11.72365705, 0.003519327,
                    12.01613674, 0.003919624,
                    12.29762575, 0.00434051,
                    12.56708522, 0.004784022,
                    12.82935694, 0.005240675,
                    13.08826821, 0.005702972,
                    13.34245457, 0.006173587,
                    13.59531574, 0.006645864,
                    13.82797848, 0.007156779,
                    14.05944231, 0.00766911,
                    14.28587854, 0.00819036,
                    14.50685953, 0.008721368,
                    14.72004211, 0.009266729,
                    14.92327099, 0.009830666,
                    15.12010698, 0.010406205,
                    15.31782083, 0.0109791,
                    15.51602442, 0.011550115,
                    15.69567496, 0.012156565,
                    15.87496664, 0.012762799,
                    16.05312135, 0.013370344,
                    16.22900547, 0.013981423,
                    16.40051069, 0.014600171,
                    16.56580098, 0.015230188,
                    16.72589258, 0.015869479,
                    16.88808634, 0.016503744,
                    17.05123224, 0.017135238,
                    17.20249998, 0.017789102,
                    17.35252354, 0.018444502,
                    17.49925873, 0.019105445,
                    17.64621938, 0.019765047,
                    17.79337514, 0.02042337,
                    17.93431888, 0.021092969,
                    18.06853085, 0.021774865,
                    18.20509134, 0.022451267,
                    18.34293087, 0.023124273,
                    18.47326374, 0.023811099,
                    18.6020398, 0.024500089,
                    18.72077992, 0.025207857,
                    18.84116452, 0.02591152,
                    18.96656381, 0.026604475,
                    19.08645174, 0.027307349,
                    19.19861785, 0.028024474,
                    19.31056919, 0.028741144,
                    19.42244357, 0.029457094,
                    19.53423905, 0.030172319,
                    19.64614387, 0.030886459,
                    19.74492998, 0.031625432,
                    19.84505806, 0.032360908,
                    19.95232331, 0.033081531,
                    20.05452074, 0.033811219,
                    20.14746825, 0.034558168,
                    20.24149331, 0.035302142,
                    20.33644216, 0.036043446,
                    20.43196586, 0.036782764,
                    20.5277897, 0.037520636,
                    20.61281675, 0.038278808,
                    20.69769345, 0.03903642,
                    20.78823157, 0.039782086,
                    20.87592273, 0.040532479,
                    20.95664499, 0.041295678,
                    21.03513484, 0.042062402,
                    21.11122572, 0.04283298,
                    21.19217463, 0.043593195,
                    21.27556842, 0.044347781,
                    21.34687384, 0.045125205,
                    21.41645983, 0.045905158,
                    21.49013154, 0.046676267,
                    21.56352165, 0.047447091,
                    21.63576365, 0.048219328,
                    21.7036602, 0.048999245,
                    21.76554835, 0.049790102,
                    21.83342207, 0.050568401,
                    21.9051746, 0.051338263,
                    21.96682138, 0.052127116,
                    22.02594212, 0.05292008,
                    22.0890376, 0.053704431,
                    22.15229605, 0.054487639,
                    22.21574697, 0.055269647,
                    22.27537873, 0.056058318,
                    22.32814311, 0.056859626,
                    22.3857899, 0.05765055,
                    22.44749408, 0.058432707,
                    22.50135065, 0.059229426,
                    22.55241117, 0.06003081,
                    22.60407869, 0.060830194,
                    22.65648813, 0.061627313,
                    22.71410964, 0.06241341,
                    22.76761149, 0.063206772,
                    22.81130477, 0.064018547,
                    22.85976873, 0.06482017,
                    22.91329469, 0.065611071,
                    22.96162214, 0.066411357,
                    23.00744641, 0.067215747,
                    23.04848121, 0.068028722,
                    23.08920017, 0.068841519,
                    23.13835765, 0.069636986,
                    23.18474079, 0.070437094,
                    23.22202816, 0.071254231,
                    23.26272496, 0.072063896,
                    23.30805509, 0.072863692,
                    23.35107702, 0.073667222,
                    23.39264311, 0.074472816,
                    23.42888167, 0.075288064,
                    23.46379829, 0.076105117,
                    23.50801996, 0.076903154,
                    23.55050312, 0.077703816,
                    23.58400368, 0.078521297,
                    23.61787894, 0.079337266,
                    23.6523985, 0.080151193,
                    23.68947425, 0.080959336,
                    23.72865459, 0.08176258,
                    23.75892805, 0.082582503,
                    23.78603085, 0.083407865,
                    23.82298486, 0.084213154,
                    23.85984574, 0.085017867,
                    23.88947382, 0.085835972,
                    23.9186609, 0.086654175,
                    23.94679817, 0.087473669,
                    23.9785529, 0.088285311,
                    24.01411358, 0.089088732,
                    24.04092141, 0.089908542,
                    24.06351936, 0.09073584,
                    24.09634303, 0.091542344,
                    24.13035312, 0.092345765,
                    24.15585461, 0.093165101,
                    24.18099108, 0.093984397,
                    24.20494282, 0.094805262,
                    24.23276068, 0.095617798,
                    24.26579565, 0.096419361,
                    24.29100718, 0.097235503,
                    24.31128426, 0.098060566,
                    24.33722303, 0.098873788,
                    24.36457535, 0.099683494,
                    24.39011651, 0.100496004,
                    24.41441344, 0.101310209,
                    24.43296897, 0.102134921,
                    24.45453935, 0.102952985,
                    24.48138454, 0.103759975,
                    24.50433827, 0.104573851,
                    24.52412476, 0.105393195,
                    24.54423044, 0.106211179,
                    24.5644478, 0.10702821,
                    24.58853274, 0.107836929,
                    24.61182847, 0.108646462,
                    24.62664948, 0.10947187,
                    24.64496104, 0.11028971,
                    24.67135029, 0.111090993,
                    24.69332602, 0.111900198,
                    24.71069756, 0.112717698,
                    24.72770021, 0.113535197,
                    24.74452339, 0.114352324,
                    24.76828232, 0.115155139,
                    24.79244306, 0.115956446,
                    24.80436952, 0.116781005,
                    24.81824263, 0.117601041,
                    24.83830019, 0.118408233,
                    24.85803271, 0.119215345,
                    24.87733092, 0.120022594,
                    24.89188861, 0.120838418,
                    24.90347834, 0.121659345,
                    24.92398239, 0.122462094,
                    24.94667088, 0.123259853,
                    24.95824878, 0.124078673,
                    24.97034122, 0.124895778,
                    24.98531204, 0.125706536,
                    25.00168119, 0.126513849,
                    25.02047973, 0.127315698,
                    25.03256594, 0.128129997,
                    25.03923363, 0.128954211,
                    25.05767426, 0.129754656,
                    25.08022208, 0.130546355,
                    25.09072635, 0.131360953,
                    25.10064504, 0.132175993,
                    25.11317815, 0.132985223,
                    25.12750263, 0.133790269,
                    25.14593951, 0.134586544,
                    25.1583689, 0.135393897,
                    25.16458718, 0.13621273,
                    25.17986593, 0.137013119,
                    25.19941592, 0.13780445,
                    25.21075787, 0.138611177,
                    25.22067676, 0.139420005,
                    25.22794995, 0.140233331,
                    25.23765744, 0.141041204,
                    25.25501397, 0.141833406,
                    25.26781331, 0.142633855,
                    25.27457227, 0.143445459,
                    25.28576614, 0.144247693,
                    25.29970103, 0.145043878,
                    25.31050753, 0.145845515,
                    25.32054535, 0.146647981,
                    25.32836349, 0.147454122,
                    25.33770305, 0.148256607,
                    25.35407718, 0.149044636,
                    25.36697556, 0.149838802,
                    25.37391555, 0.150643973,
                    25.38260241, 0.151445051,
                    25.39268221, 0.152242732,
                    25.40523706, 0.153034894,
                    25.4186373, 0.153824734,
                    25.42377102, 0.154630106,
                    25.42962453, 0.155433403,
                    25.44545698, 0.156216484,
                    25.45902415, 0.157003341,
                    25.46747551, 0.157799562,
                    25.47461536, 0.158597691,
                    25.48041155, 0.159397795,
                    25.49146608, 0.160186938,
                    25.50496917, 0.160970625,
                    25.50950115, 0.161771235,
                    25.51317675, 0.162572869,
                    25.5252606, 0.163357373,
                    25.53621704, 0.164143434,
                    25.54372326, 0.164935604,
                    25.55054801, 0.165728458,
                    25.55648386, 0.166522405,
                    25.56862013, 0.167303554,
                    25.58455631, 0.168076611,
                    25.58902646, 0.168871487,
                    25.59074401, 0.169671112,
                    25.60172838, 0.170451934,
                    25.61249844, 0.171232534,
                    25.62217876, 0.172014627,
                    25.62915095, 0.172801384,
                    25.6315321, 0.173596497,
                    25.63960766, 0.174379814,
                    25.65219237, 0.175153657,
                    25.65953683, 0.175937131,
                    25.66507978, 0.176723498,
                    25.66951969, 0.177511392,
                    25.67495333, 0.178296705,
                    25.68830415, 0.179065871,
                    25.69742512, 0.179842694,
                    25.69708638, 0.180637421,
                    25.70334941, 0.181418582,
                    25.71632018, 0.182185971,
                    25.72278429, 0.182965482,
                    25.72622706, 0.183750285,
                    25.73187946, 0.184530133,
                    25.73809904, 0.185308244,
                    25.7484318, 0.186077672,
                    25.75645313, 0.186851005,
                    25.7574575, 0.187637466,
                    25.76383961, 0.188412769,
                    25.77723526, 0.189173708,
                    25.78514473, 0.189944779,
                    25.78971093, 0.190721782,
                    25.79437652, 0.191497973,
                    25.79904673, 0.192273538,
                    25.80910751, 0.193037926,
                    25.81774629, 0.193804486,
                    25.82012772, 0.194582694,
                    25.82374224, 0.195357872,
                    25.82945227, 0.196128334,
                    25.83566525, 0.196897199,
                    25.84226865, 0.19766469,
                    25.84385485, 0.198441422,
                    25.84372995, 0.19922088,
                    25.85352605, 0.199980292,
                    25.863009, 0.200739712,
                    25.86552518, 0.201512177,
                    25.86866223, 0.202282823,
                    25.87320888, 0.203050103,
                    25.88064136, 0.203811127,
                    25.89098013, 0.204565856,
                    25.89399698, 0.20533433,
                    25.89364752, 0.206108802,
                    25.90368089, 0.206862331,
                    25.91466903, 0.207613392,
                    25.91756744, 0.208379708,
                    25.92011498, 0.209146115,
                    25.92165823, 0.209913896,]
        pp_table = np.array(pp_table).reshape(-1, 2).T.tolist()
        pe_table = [6.201145019, 0,
                    7.406784123, 0.000222709,
                    7.785127381, 0.000466271,
                    8.138957899, 0.000756893,
                    8.512860278, 0.001007219,
                    8.923006579, 0.001185563,
                    9.30323792, 0.001421559,
                    9.658119036, 0.001706266,
                    10.03280602, 0.001951205,
                    10.41534024, 0.00217981,
                    10.76908969, 0.002463862,
                    11.11595909, 0.002760439,
                    11.42160461, 0.003136838,
                    11.72365705, 0.003519327,
                    12.01613674, 0.003919624,
                    12.29762575, 0.00434051,
                    12.56708522, 0.004784022,
                    12.82935694, 0.005240675,
                    13.08826821, 0.005702972,
                    13.34245457, 0.006173587,
                    13.59531574, 0.006645864,
                    13.82797848, 0.007156779,
                    14.05944231, 0.00766911,
                    14.28587854, 0.00819036,
                    14.50685953, 0.008721368,
                    14.72004211, 0.009266729,
                    14.92327099, 0.009830666,
                    15.12010698, 0.010406205,
                    15.31782083, 0.0109791,
                    15.51602442, 0.011550115,
                    15.69567496, 0.012156565,
                    15.87496664, 0.012762799,
                    16.05312135, 0.013370344,
                    16.22900547, 0.013981423,
                    16.40051069, 0.014600171,
                    16.56580098, 0.015230188,
                    16.72589258, 0.015869479,
                    16.88808634, 0.016503744,
                    17.05123224, 0.017135238,
                    17.20249998, 0.017789102,
                    17.35252354, 0.018444502,
                    17.49925873, 0.019105445,
                    17.64621938, 0.019765047,
                    17.79337514, 0.02042337,
                    17.93431888, 0.021092969,
                    18.06853085, 0.021774865,
                    18.20509134, 0.022451267,
                    18.34293087, 0.023124273,
                    18.47326374, 0.023811099,
                    18.6020398, 0.024500089,
                    18.72077992, 0.025207857,
                    18.84116452, 0.02591152,
                    18.96656381, 0.026604475,
                    19.08645174, 0.027307349,
                    19.19861785, 0.028024474,
                    19.31056919, 0.028741144,
                    19.42244357, 0.029457094,
                    19.53423905, 0.030172319,
                    19.64614387, 0.030886459,
                    19.74492998, 0.031625432,
                    19.84505806, 0.032360908,
                    19.95232331, 0.033081531,
                    20.05452074, 0.033811219,
                    20.14746825, 0.034558168,
                    20.24149331, 0.035302142,
                    20.33644216, 0.036043446,
                    20.43196586, 0.036782764,
                    20.5277897, 0.037520636,
                    20.61281675, 0.038278808,
                    20.69769345, 0.03903642,
                    20.78823157, 0.039782086,
                    20.87592273, 0.040532479,
                    20.95664499, 0.041295678,
                    21.03513484, 0.042062402,
                    21.11122572, 0.04283298,
                    21.19217463, 0.043593195,
                    21.27556842, 0.044347781,
                    21.34687384, 0.045125205,
                    21.41645983, 0.045905158,
                    21.49013154, 0.046676267,
                    21.56352165, 0.047447091,
                    21.63576365, 0.048219328,
                    21.7036602, 0.048999245,
                    21.76554835, 0.049790102,
                    21.83342207, 0.050568401,
                    21.9051746, 0.051338263,
                    21.96682138, 0.052127116,
                    22.02594212, 0.05292008,
                    22.0890376, 0.053704431,
                    22.15229605, 0.054487639,
                    22.21574697, 0.055269647,
                    22.27537873, 0.056058318,
                    22.32814311, 0.056859626,
                    22.3857899, 0.05765055,
                    22.44749408, 0.058432707,
                    22.50135065, 0.059229426,
                    22.55241117, 0.06003081,
                    22.60407869, 0.060830194,
                    22.65648813, 0.061627313,
                    22.71410964, 0.06241341,
                    22.76761149, 0.063206772,
                    22.81130477, 0.064018547,
                    22.85976873, 0.06482017,
                    22.91329469, 0.065611071,
                    22.96162214, 0.066411357,
                    23.00744641, 0.067215747,
                    23.04848121, 0.068028722,
                    23.08920017, 0.068841519,
                    23.13835765, 0.069636986,
                    23.18474079, 0.070437094,
                    23.22202816, 0.071254231,
                    23.26272496, 0.072063896,
                    23.30805509, 0.072863692,
                    23.35107702, 0.073667222,
                    23.39264311, 0.074472816,
                    23.42888167, 0.075288064,
                    23.46379829, 0.076105117,
                    23.50801996, 0.076903154,
                    23.55050312, 0.077703816,
                    23.58400368, 0.078521297,
                    23.61787894, 0.079337266,
                    23.6523985, 0.080151193,
                    23.68947425, 0.080959336,
                    23.72865459, 0.08176258,
                    23.75892805, 0.082582503,
                    23.78603085, 0.083407865,
                    23.82298486, 0.084213154,
                    23.85984574, 0.085017867,
                    23.88947382, 0.085835972,
                    23.9186609, 0.086654175,
                    23.94679817, 0.087473669,
                    23.9785529, 0.088285311,
                    24.01411358, 0.089088732,
                    24.04092141, 0.089908542,
                    24.06351936, 0.09073584,
                    24.09634303, 0.091542344,
                    24.13035312, 0.092345765,
                    24.15585461, 0.093165101,
                    24.18099108, 0.093984397,
                    24.20494282, 0.094805262,
                    24.23276068, 0.095617798,
                    24.26579565, 0.096419361,
                    24.29100718, 0.097235503,
                    24.31128426, 0.098060566,
                    24.33722303, 0.098873788,
                    24.36457535, 0.099683494,
                    24.39011651, 0.100496004,
                    24.41441344, 0.101310209,
                    24.43296897, 0.102134921,
                    24.45453935, 0.102952985,
                    24.48138454, 0.103759975,
                    24.50433827, 0.104573851,
                    24.52412476, 0.105393195,
                    24.54423044, 0.106211179,
                    24.5644478, 0.10702821,
                    24.58853274, 0.107836929,
                    24.61182847, 0.108646462,
                    24.62664948, 0.10947187,
                    24.64496104, 0.11028971,
                    24.67135029, 0.111090993,
                    24.69332602, 0.111900198,
                    24.71069756, 0.112717698,
                    24.72770021, 0.113535197,
                    24.74452339, 0.114352324,
                    24.76828232, 0.115155139,
                    24.79244306, 0.115956446,
                    24.80436952, 0.116781005,
                    24.81824263, 0.117601041,
                    24.83830019, 0.118408233,
                    24.85803271, 0.119215345,
                    24.87733092, 0.120022594,
                    24.89188861, 0.120838418,
                    24.90347834, 0.121659345,
                    24.92398239, 0.122462094,
                    24.94667088, 0.123259853,
                    24.95824878, 0.124078673,
                    24.97034122, 0.124895778,
                    24.98531204, 0.125706536,
                    25.00168119, 0.126513849,
                    25.02047973, 0.127315698,
                    25.03256594, 0.128129997,
                    25.03923363, 0.128954211,
                    25.05767426, 0.129754656,
                    25.08022208, 0.130546355,
                    25.09072635, 0.131360953,
                    25.10064504, 0.132175993,
                    25.11317815, 0.132985223,
                    25.12750263, 0.133790269,
                    25.14593951, 0.134586544,
                    25.1583689, 0.135393897,
                    25.16458718, 0.13621273,
                    25.17986593, 0.137013119,
                    25.19941592, 0.13780445,
                    25.21075787, 0.138611177,
                    25.22067676, 0.139420005,
                    25.22794995, 0.140233331,
                    25.23765744, 0.141041204,
                    25.25501397, 0.141833406,
                    25.26781331, 0.142633855,
                    25.27457227, 0.143445459,
                    25.28576614, 0.144247693,
                    25.29970103, 0.145043878,
                    25.31050753, 0.145845515,
                    25.32054535, 0.146647981,
                    25.32836349, 0.147454122,
                    25.33770305, 0.148256607,
                    25.35407718, 0.149044636,
                    25.36697556, 0.149838802,
                    25.37391555, 0.150643973,
                    25.38260241, 0.151445051,
                    25.39268221, 0.152242732,
                    25.40523706, 0.153034894,
                    25.4186373, 0.153824734,
                    25.42377102, 0.154630106,
                    25.42962453, 0.155433403,
                    25.44545698, 0.156216484,
                    25.45902415, 0.157003341,
                    25.46747551, 0.157799562,
                    25.47461536, 0.158597691,
                    25.48041155, 0.159397795,
                    25.49146608, 0.160186938,
                    25.50496917, 0.160970625,
                    25.50950115, 0.161771235,
                    25.51317675, 0.162572869,
                    25.5252606, 0.163357373,
                    25.53621704, 0.164143434,
                    25.54372326, 0.164935604,
                    25.55054801, 0.165728458,
                    25.55648386, 0.166522405,
                    25.56862013, 0.167303554,
                    25.58455631, 0.168076611,
                    25.58902646, 0.168871487,
                    25.59074401, 0.169671112,
                    25.60172838, 0.170451934,
                    25.61249844, 0.171232534,
                    25.62217876, 0.172014627,
                    25.62915095, 0.172801384,
                    25.6315321, 0.173596497,
                    25.63960766, 0.174379814,
                    25.65219237, 0.175153657,
                    25.65953683, 0.175937131,
                    25.66507978, 0.176723498,
                    25.66951969, 0.177511392,
                    25.67495333, 0.178296705,
                    25.68830415, 0.179065871,
                    25.69742512, 0.179842694,
                    25.69708638, 0.180637421,
                    25.70334941, 0.181418582,
                    25.71632018, 0.182185971,
                    25.72278429, 0.182965482,
                    25.72622706, 0.183750285,
                    25.73187946, 0.184530133,
                    25.73809904, 0.185308244,
                    25.7484318, 0.186077672,
                    25.75645313, 0.186851005,
                    25.7574575, 0.187637466,
                    25.76383961, 0.188412769,
                    25.77723526, 0.189173708,
                    25.78514473, 0.189944779,
                    25.78971093, 0.190721782,
                    25.79437652, 0.191497973,
                    25.79904673, 0.192273538,
                    25.80910751, 0.193037926,
                    25.81774629, 0.193804486,
                    25.82012772, 0.194582694,
                    25.82374224, 0.195357872,
                    25.82945227, 0.196128334,
                    25.83566525, 0.196897199,
                    25.84226865, 0.19766469,
                    25.84385485, 0.198441422,
                    25.84372995, 0.19922088,
                    25.85352605, 0.199980292,
                    25.863009, 0.200739712,
                    25.86552518, 0.201512177,
                    25.86866223, 0.202282823,
                    25.87320888, 0.203050103,
                    25.88064136, 0.203811127,
                    25.89098013, 0.204565856,
                    25.89399698, 0.20533433,
                    25.89364752, 0.206108802,
                    25.90368089, 0.206862331,
                    25.91466903, 0.207613392,
                    25.91756744, 0.208379708,
                    25.92011498, 0.209146115,
                    25.921658, 0.209913896,]
        pe_table = np.array(pe_table).reshape(-1, 2).T.tolist()

        # create hardening table for matrix
        self.hardening_table_matrix = pe_table
        # create hardening table for fiber
        self.hardening_table_fiber = pp_table

        self.sim_paras = {
            "size": size,
            "radius_mu": radius_mu,
            "radius_std": radius_std,
            "vol_req": vol_req,
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
            "num_cpu": num_cpu}

        self.mini_dist_factor = min_dist_factor

        # print simulation information to screen
        if print_info:
            self._print_sim_info(info=self.sim_paras)

    def _get_sim_info(self) -> None:
        """get simulation information"""
        self.sim_info = {
            "job_name": "J2StrucMesh",
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
            "hardening_table_fiber": self.hardening_table_fiber,
            "hardening_table_matrix": self.hardening_table_matrix,
            "num_steps": self.num_steps,
            "simulation_time": self.simulation_time,
            "strain": self.strain,
            "num_cpu": self.num_cpu,
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
        # create microstructure
        self.microstructure = CircleParticles(
            length=self.size,
            width=self.size,
            radius_mu=self.radius_mu,
            radius_std=self.radius_std,
            vol_req=self.vol_req,
            dist_min_factor=self.mini_dist_factor,
        )
        self.microstructure.generate_microstructure(seed=self.seed)
        self.microstructure.to_abaqus_format()
        self.microstructure.plot_microstructure(save_figure=True,
                                                fig_name="rve_{}.png".
                                                format(self.seed))
        # save the microstructure to rgmsh.npy format
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
