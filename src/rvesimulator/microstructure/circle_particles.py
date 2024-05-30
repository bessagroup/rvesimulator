"""microstructure generator for 2D RVE with different size of disks/circles"""
#                                                                       Modules
# =============================================================================
# standard
import json
import logging
import math
import time
from typing import Any, Tuple

# Third party
import numpy as np
from scipy.spatial import distance_matrix

# local
from .base import MicrostructureGenerator

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class CircleParticles(MicrostructureGenerator):
    """2D RVE with different size of disks/circles

    Parameters
    ----------
    MicrostructureGenerator : class
        parent class of microstructure generator

    Examples
    --------
    >>> from rvesimulator.microstructure import CircleParticles
    >>> circle_particles = CircleParticles(
    ...     length=1,
    ...     width=1,
    ...     radius_mu=0.1,
    ...     radius_std=0.01,
    ...     vol_req=0.1,
    ... )
    >>> circle_particles.generate_microstructure()
    >>> circle_particles.plot_microstructure()
    >>> circle_particles.to_abaqus_format()
    >>> circle_particles.crate_rgmsh()
    """

    def __init__(
        self,
        length: float,
        width: float,
        radius_mu: float,
        radius_std: float,
        vol_req: float,
        num_guess_max: int = 1000000,
        num_inclusion_max: int = 750,
        num_cycle_max: int = 15,
        dist_min_factor: float = 1.1,
        stirring_iters: int = 100,
        print_log: bool = False,
        vol_frac_tol: float = 0.01,
        distance_tol: float = 0.02,
    ) -> None:
        """Initialization

        Parameters
        ----------
        length : float
            length of RVE
        width : float
            width of RVE
        radius_mu : float
            mean of circle's radius
        radius_std : float
            std of circle's radius
        vol_req : float
            required volume fraction
        num_guess_max : int, optional
            maximum guess for inclusions, by default 50000
        num_inclusion_max : int, optional
            maximum inclusion inside RVE, by default 750
        num_cycle_max : int, optional
            iteration cycles, by default 15
        dist_min_factor : float, optional
            distance factor, by default 2.07
        """
        # geometry information of the 2D RVE with homogeneous circles
        self.length = length
        self.width = width
        self.radius_mu = radius_mu
        self.radius_std = radius_std
        self.vol_req = vol_req
        self.print_log = print_log

        # Initialization of the algorithm
        self.dist_min_factor = dist_min_factor
        self.num_guess_max = num_guess_max
        self.num_inclusions_max = num_inclusion_max
        self.num_cycles_max = num_cycle_max
        self.stirring_iters = stirring_iters
        self.vol_frac_tol = vol_frac_tol
        self.distance_tol = distance_tol

        # adjust the radius_mu for achieving the required volume fraction, only when the radius_std is zero:
        if self.radius_std == 0:
            num_particles = np.round(vol_req*length*width/(np.pi*radius_mu**2))
            self.min_dis_threshold = max(np.sqrt(length*width/num_particles) - radius_mu, 2*radius_mu)
            current_vfrac = num_particles*np.pi*radius_mu**2/(length*width)
            if abs(current_vfrac-vol_req) > vol_frac_tol:
                radius_mu = np.sqrt(vol_req * length * width / (num_particles * np.pi))
            self.radius_mu = radius_mu
        else:
            num_particles = np.round(vol_req*length*width/(np.pi*radius_mu**2))
            self.min_dis_threshold = None
            self.dist_min_factor = 1 + 2*0.01*(np.sqrt(length*width/num_particles) - 0.5*radius_mu)/(radius_mu*vol_req)
            print(f"dist_min_factor: {self.dist_min_factor}")



    def _parameter_initialization(self) -> None:
        """Initialize the parameters"""
        self.num_change = 3
        self.num_cycle = 0
        self.vol_frac = 0
        self.vol_total = self.length * self.width

        # initial coordinate for position of inclusion
        self.len_start = -1 * self.radius_mu
        self.len_end = self.length + self.radius_mu
        self.wid_start = -1 * self.radius_mu
        self.wid_end = self.width + self.radius_mu

        # inclusion location is a nx4 numpy array x, y, r, p (partition)
        self.inclusion_positions = None

    def generate_microstructure(
        self,
        seed: Any = None,
    ) -> None:
        """generate microstructure

        Parameters
        ----------
        seed : Any, optional
            seed number, by default None
        """

        # decide to use seed or not

        self.rng = np.random.default_rng(seed=seed)
        # counting time generating an RVE
        logging.basicConfig(level=logging.INFO,
                            filename='rve_simulation.log', filemode='w')
        self.logger = logging.getLogger("microstructure")
        # Create a buffer handler
        self.logger.info("==================================================")
        self.logger.info("Start generating microstructure")
        self.logger.info(f"seed: {seed}")
        self.logger.info(f"length: {self.length}")
        self.logger.info(f"width: {self.width}")
        self.logger.info(f"radius_mu: {self.radius_mu}")
        self.logger.info(f"radius_std: {self.radius_std}")
        self.logger.info(f"vol_req: {self.vol_req}")
        start_time = time.time()
        self._parameter_initialization()
        self._procedure_initialization()
        self._core_iteration()
        end_time = time.time()
        self.time_usage = end_time - start_time
        self.logger.info(f"time usage: {self.time_usage}")
        self.logger.info("End generating microstructure")
        self.logger.info("==================================================")
        if self.print_log:
            file_handler = logging.FileHandler("microstructure.log")
            file_handler.setLevel(logging.INFO)
            # Create a logging format
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            file_handler.setFormatter(formatter)
            # Add the handlers to the logger
            self.logger.addHandler(file_handler)

        # get microstructure info
        self.microstructure_info = {
            "inclusion_location_information": self.inclusion_positions.tolist(),
            "radius_mu": self.radius_mu,
            "radius_std": self.radius_std,
            "length_start": self.len_start,
            "width_start": self.wid_start,
            "length_end": self.len_end,
            "width_end": self.wid_end,
            "vol_frac": self.vol_frac,
            "seed": seed,
        }


    def to_abaqus_format(
        self, file_name: str = "micro_structure_info.json"
    ) -> dict:
        """save rve microstructure to abaqus format

        Parameters
        ----------
        file_name : str, optional
            name of saved file, by default "micro_structure_info.json"

        Returns
        -------
        dict
            a dict contains all info of rve microstructure
        """
        with open(file_name, "w") as fp:
            json.dump(self.microstructure_info, fp)


    def plot_microstructure(
        self,
        save_figure: bool = False,
        fig_name: str = "mircostructure.png",
        **kwargs,
    ) -> None:
        """plot microstructure

        Parameters
        ----------
        save_figure : bool, optional
            save figure or not, by default False
        fig_name : str, optional
            name of figure, by default "mircostructure.png"
        """
        self.circle_plot(
            inclusions=self.inclusion_positions,
            length=self.length,
            width=self.width,
            vol_frac=self.vol_frac,
            save_figure=save_figure,
            fig_name=fig_name,
            **kwargs,
        )

    def _procedure_initialization(self) -> None:
        """This function is used to generate the first disk and
        assign the initial values of the algorithm
        """

        # initialization(generate the first inclusion randomly)
        self.inclusion_min_dis_vector = np.zeros(
            (self.num_inclusions_max, self.num_cycles_max + 1, 2)
        )
        self.num_inclusions = 1
        # generate the location of the first inclusion
        # the first inclusion is generated with one partition
        inclusion_temp = self.generate_random_inclusions(
            len_start=self.radius_mu + 2*self.distance_tol,
            len_end=self.length - self.radius_mu - 2*self.distance_tol,
            wid_start=self.radius_mu + 2*self.distance_tol,
            wid_end=self.width - self.radius_mu - 2*self.distance_tol,
            radius_mu=self.radius_mu,
            radius_std=0,
            rng=self.rng,
        )
        # update the volume fraction information
        self.vol_frac = self.inclusion_volume(self.radius_mu) / self.vol_total
        self.inclusion_positions = np.zeros((1, 4))
        self.inclusion_positions[0, 0:3] = inclusion_temp.T
        self.inclusion_positions[0, 3] = 1

    def _core_iteration(self) -> None:
        """core iteration part of the micro-structure generation method"""

        # the main loop of the algorithm
        while (
            abs(self.vol_frac - self.vol_req) > self.vol_frac_tol
            and self.num_cycle < self.num_cycles_max
        ):
            # ================================================================#
            #                   generate the inclusions randomly                  #
            # ================================================================#
            self.num_trial = 1
            while (
                self.num_trial < self.num_guess_max
                and abs(self.vol_frac - self.vol_req) > self.vol_frac_tol
                and self.num_inclusions < self.num_inclusions_max
            ):
                # update the info of number trial
                self.num_trial = self.num_trial + 1
                inclusion_temp = self.generate_random_inclusions(
                    len_start=0,
                    len_end=self.length,
                    wid_start=0,
                    wid_end=self.width,
                    radius_mu=self.radius_mu,
                    radius_std=self.radius_std,
                    rng=self.rng,
                )
                # check the location of the inclusion and
                new_inclusion = self.new_positions(
                    x_center=inclusion_temp[0, 0],
                    y_center=inclusion_temp[1, 0],
                    radius=inclusion_temp[2, 0],
                    length=self.length,
                    width=self.width,
                )
                if new_inclusion[0, 3] == 4:
                    self.logger.info("generate inclusion, vertex check ...")
                    # if the temp inclusion locates at un-proper location for mesh
                    while self.vertices_mesh_loc(new_inclusion) == "fail":
                        inclusion_temp = self.generate_random_inclusions(
                            len_start=0,
                            len_end=self.length,
                            wid_start=0,
                            wid_end=self.width,
                            radius_mu=self.radius_mu,
                            radius_std=self.radius_std,
                            rng=self.rng,
                        )
                        new_inclusion = self.new_positions(
                            x_center=inclusion_temp[0, 0],
                            y_center=inclusion_temp[1, 0],
                            radius=inclusion_temp[2, 0],
                            length=self.length,
                            width=self.width,
                        )
                    self.logger.info("generate inclusion, vertex check pass")
                elif new_inclusion[0, 3] == 2 or new_inclusion[0, 3] == 1:
                    self.logger.info("generate inclusion, edge check ...")
                    # if the temp inclusion locates at un-proper location for mesh
                    while self.proper_edge_mesh_location(new_inclusion) == "fail":
                        inclusion_temp = self.generate_random_inclusions(
                            len_start=0,
                            len_end=self.length,
                            wid_start=0,
                            wid_end=self.width,
                            radius_mu=self.radius_mu,
                            radius_std=self.radius_std,
                            rng=self.rng,
                        )
                        new_inclusion = self.new_positions(
                            x_center=inclusion_temp[0, 0],
                            y_center=inclusion_temp[1, 0],
                            radius=inclusion_temp[2, 0],
                            length=self.length,
                            width=self.width,
                        )
                    self.logger.info("generate inclusion, edge check pass")
                # check the overlap of new inclusion
                overlap_status = self.overlap_check(
                    new_inclusion=new_inclusion,
                    inclusion_pos=self.inclusion_positions.copy(),
                    dist_factor=self.dist_min_factor,
                    distance_tol=self.distance_tol,
                    min_dis_threshold=self.min_dis_threshold
                )
                if overlap_status == 0:
                    self.inclusion_positions = np.vstack(
                        (self.inclusion_positions, new_inclusion)
                    )
                    self.vol_frac = (
                        self.vol_frac
                        + self.inclusion_volume(new_inclusion[0, 2]) / self.vol_total
                    )
                    self.num_inclusions = self.num_inclusions + new_inclusion.shape[0]
                del new_inclusion

            # ================================================================#
            #                   stirring the inclusions (first stage)             #
            # ================================================================#
            ii = 0
            if self.inclusion_positions.shape[0] < self.num_inclusions_max:
                # for every point, stirring is needed!
                while ii < self.inclusion_positions.shape[0]:
                    (
                        self.inclusion_min_dis_vector,
                        min_index,
                        min_dis,
                    ) = self.min_dis_index(
                        self.inclusion_positions[ii, 0:2],
                        self.inclusion_positions.copy(),
                        self.inclusion_min_dis_vector,
                        ii,
                        self.num_cycle,
                    )
                    # generate the new inclusion location
                    new_inclusion_temp = self.gen_heuristic_inclusions(
                        ref_point=self.inclusion_positions[min_index, 0:3].copy(),
                        inclusion_temp=self.inclusion_positions[ii, 0:3].copy(),
                        dist_factor=self.dist_min_factor,
                        rng=self.rng,
                    )
                    # check the overlap of new inclusion
                    new_inclusion = self.new_positions(
                        x_center=new_inclusion_temp[0, 0],
                        y_center=new_inclusion_temp[0, 1],
                        radius=new_inclusion_temp[0, 2],
                        length=self.length,
                        width=self.width,
                    )
                    # check proper location for mesh
                    # max stirring iteration
                    stirring_iter = 0
                    if new_inclusion[0, 3] == 4:
                        self.logger.info("stirring inclusion,vertex check ...")
                        while (
                            self.vertices_mesh_loc(new_inclusion) == "fail"
                            and stirring_iter < self.stirring_iters
                        ):
                            # generate new inclusion
                            self.logger.info(f"iter: {stirring_iter}")
                            new_inclusion_temp = self.gen_heuristic_inclusions(
                                ref_point=self.inclusion_positions[
                                    min_index, 0:3
                                ].copy(),
                                inclusion_temp=self.inclusion_positions[
                                    ii, 0:3].copy(),
                                dist_factor=self.dist_min_factor,
                                rng=self.rng,
                            )
                            # check the overlap of new inclusion
                            new_inclusion = self.new_positions(
                                x_center=new_inclusion_temp[0, 0],
                                y_center=new_inclusion_temp[0, 1],
                                radius=new_inclusion_temp[0, 2],
                                length=self.length,
                                width=self.width,
                            )
                            stirring_iter += 1
                        if stirring_iter == self.stirring_iters:
                            self.logger.error(
                                "stirring vertex check failed")
                        else:
                            self.logger.info("stirring vertex check pass")
                    elif new_inclusion[0, 3] == 2:
                        self.logger.info("stirring inclusion, edge check ...")
                        # check proper location for mesh
                        while self.proper_edge_mesh_location(
                            new_inclusion
                        ) == "fail" and stirring_iter < self.stirring_iters:
                            # logger
                            self.logger.info(f"iter: {stirring_iter}")
                            new_inclusion_temp = self.gen_heuristic_inclusions(
                                ref_point=self.inclusion_positions[
                                    min_index, 0:3
                                ].copy(),
                                inclusion_temp=self.inclusion_positions[
                                    ii, 0:3].copy(),
                                dist_factor=self.dist_min_factor,
                                rng=self.rng,
                            )
                            # check the overlap of new inclusion
                            new_inclusion = self.new_positions(
                                x_center=new_inclusion_temp[0, 0],
                                y_center=new_inclusion_temp[0, 1],
                                radius=new_inclusion_temp[0, 2],
                                length=self.length,
                                width=self.width,
                            )
                            stirring_iter += 1
                        if stirring_iter == self.stirring_iters:
                            self.logger.error("stirring edge check failed")
                        # else:
                        #     self.logger.info("stirring edge check pass")
                    overlap_status = self.overlap_check(
                        new_inclusion=new_inclusion,
                        inclusion_pos=self.inclusion_positions.copy(),
                        dist_factor=self.dist_min_factor,
                        stage="step_two",
                        inclusion_index=ii,
                        distance_tol=self.distance_tol,
                        min_dis_threshold=self.min_dis_threshold,
                    )
                    # check: if the new inclusions(cause it maybe more than
                    # 1 inclusion centers) will overlap with the
                    # remaining ones or not
                    if overlap_status == 0:
                        ii = self._update_inclusion_position(
                            new_inclusion=new_inclusion, iter=ii
                        )
                    else:
                        ii = ii + int(self.inclusion_positions[ii, 3])

                    del new_inclusion, new_inclusion_temp
            # end of one cycle
            self.num_cycle = self.num_cycle + 1
            if self.num_cycle == self.num_cycles_max:
                self.logger.error("The algorithm reached the maximum cycle number")

    def _update_inclusion_position(self, new_inclusion: np.ndarray, iter: int) -> int:
        """update the inclusion position

        Parameters
        ----------
        new_inclusion : np.ndarray
            the generated new inclusion
        iter : int
            determine the nex inclusion should be analysis

        Returns
        -------
        iter : int
            he updated number of index
        """
        # check the location compatibility
        if new_inclusion[0, 3] != new_inclusion.shape[0]:
            raise ValueError("inclusion number compatibility issue")
        inclusion_portion = int(self.inclusion_positions[iter, 3].copy())
        self.inclusion_positions = np.delete(
            self.inclusion_positions,
            tuple(i + iter for i in range(inclusion_portion)),
            axis=0,
        )
        self.inclusion_positions = np.insert(
            self.inclusion_positions, (iter), new_inclusion, axis=0
        )
        iter = iter + int(new_inclusion[0, 3])
        assert type(iter) == int

        return iter

    def crate_rgmsh(self, num_discrete: int = 10) -> np.ndarray:
        """create rgmsh numpy array for crate

        Parameters
        ----------
        num_discrete : int, optional
            number of discrete partition, by default 10

        Returns
        -------
        np.ndarray
            2d numpy array that contains the micro-structure information
        """

        self.rgmsh = np.zeros((num_discrete, num_discrete))
        grid_len = self.length / num_discrete
        grid_wid = self.width / num_discrete
        radius = self.inclusion_positions[:, 2].reshape(-1, 1)
        for ii in range(num_discrete):
            for jj in range(num_discrete):
                loc_temp = np.array(
                    [
                        [
                            self.length / (2 * num_discrete) + ii * grid_len,
                            self.width / (2 * num_discrete) + jj * grid_wid,
                        ]
                    ]
                )
                # distance measure
                points_dis_temp = distance_matrix(
                    self.inclusion_positions[:, 0:2],
                    loc_temp,
                )

                if (points_dis_temp - radius).min() < 0:
                    self.rgmsh[ii, jj] = 1

        return self.rgmsh.T

    def vertices_mesh_loc(self, inclusion: np.ndarray) -> int:
        """identify proper vertices location for meshing

        Parameters
        ----------
        inclusion : np.ndarray
            temp inclusion

        Returns
        -------
        int
            status of the inclusion(0: improper, 1: proper)
        """
        # # reformat the inclusion location
        # inclusion = inclusion.reshape((-1, 4))
        # vertices = np.array([[0, 0],
        #                      [0, self.width],
        #                      [self.length, self.width],           
        #                      [self.length, 0]])
        # # calculate the distance between the inclusion and the vertices
        # points_dis_temp = distance_matrix(
        #     vertices,
        #     inclusion[:, 0:2],
        # )
        # min_points_dis = points_dis_temp.min()
        # if 0.95*inclusion[0, 2] < min_points_dis < np.sqrt(2)*inclusion[0, 2]:
        #     return "fail"
        # else:
        #     return "pass"
        return self.proper_edge_mesh_location(inclusion)

    def proper_edge_mesh_location(self, inclusion: np.ndarray) -> int:
        """identify proper edge location for meshing

        Parameters
        ----------
        inclusion : np.ndarray
            temp inclusion

        Returns
        -------
        int
            status of the inclusion(0: improper, 1: proper)
        """
        
        # reformat the inclusion location
        inclusion = inclusion.reshape((-1, 4))

        # ensure the radius of the inclusion is 2 times larger than the edge tolerance
        radius = inclusion[0, 2]
        if radius < 4*self.distance_tol:
            return "fail"
        
        status_x = "fail"
        # for x edges
        dis_x = np.abs(np.array([inclusion[:, 0], self.width - inclusion[:, 0]]))
        # check minimum distance in the inside region
        if dis_x.min() > radius + 2*self.distance_tol:
            status_x = "pass"

        # making exception for the volume fraction larger than 0.4
        # if self.vol_req > 0.4:
        # check the minimum distance in the boundary region
            # if dis_x.min() < radius/2.0:
            #     status_x = "pass"

        # if 0.95*inclusion[0, 2] < dis_x.min() < inclusion[0, 2]:
        #     return "fail"
        # elif inclusion[0, 2] < dis_x.min() < 1.05*inclusion[0, 2]:
        #     return "fail"
        # for y edges

        status_y = "fail"
        dis_y = np.abs(np.array([inclusion[:, 1], self.length - inclusion[:, 1]]))
        # check minimum distance in the inside region
        if dis_y.min() > radius + 2*self.distance_tol:
            status_y = "pass"

        # making exception for the volume fraction larger than 0.4
        # if self.vol_req > 0.4:
        # check the minimum distance in the boundary region
            # if dis_y.min() < radius/2.0:
            #     status_y = "pass"

        if status_x == "pass" and status_y == "pass":
            return "pass"
        else:
            return "fail"

        # if 0.95*inclusion[0, 2] < dis_y.min() < inclusion[0, 2]:
        #     return "fail"
        # elif inclusion[0, 2] < dis_y.min() < 1.05*inclusion[0, 2]:
        #     return 0
        # return 

    def new_positions(
        self,
        x_center: float,
        y_center: float,
        radius: float,
        width: float,
        length: float,
    ) -> np.ndarray:
        """This is the function used to determine the locations of disks
        considering the boundary of RVE. To be specific, the disk should
        be divided into two parts or four parts if it is located in the
        outskirt area

        Parameters
        ----------
        x_center : float
            new center of X axis
        y_center : float
            new center of Y axis
        radius : float
            radius of the disk
        width : float
            width of RVE
        length : float
            length of RVE

        Returns
        -------
        np.ndarray
            XC, YC, split  which is the new locations of this inclusion
            (because in some locations, the circles need to be split)
        ####################################
        #     #          2_2       #       #
        # 4_1 #                    #4_2    #
        ####################################
        #     #                    #       #
        # 2_3 #          1         #   2_4 #
        #     #                    #       #
        #     #                    #       #
        ####################################
        # 4_3 #         2_1        # 4_4   #
        #     #                    #       #
        ####################################

        """

        new_inclusion = np.zeros((1, 4))
        if (
            radius <= x_center <= length - radius
            and radius <= y_center <= width - radius
        ):
            # locate in center region and split = 1
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 1)

        elif length - radius > x_center > radius > y_center:
            # location 2_1
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 2)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center, y_center + width, radius, 2]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif radius < x_center < length - radius and y_center > width - radius:
            # location 2_2
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 2)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center, y_center - width, radius, 2]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif width - radius > y_center > radius > x_center:
            # location 2_3
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 2)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center + length, y_center, radius, 2]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif radius < y_center < width - radius and x_center > length - radius:
            # location 2_4
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 2)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center - length, y_center, radius, 2]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif x_center < radius and y_center > width - radius:
            # location 4_1
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 4)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center + length, y_center, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array(
                        [x_center + length, y_center - width, radius, 4]
                    ).reshape((1, 4)),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center, y_center - width, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif x_center > length - radius and y_center > width - radius:
            # location 4_2
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 4)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center - length, y_center, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array(
                        [x_center - length, y_center - width, radius, 4]
                    ).reshape((1, 4)),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center, y_center - width, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif x_center < radius and y_center < radius:
            # location 4_3
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 4)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center + length, y_center, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array(
                        [x_center + length, y_center + width, radius, 4]
                    ).reshape((1, 4)),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center, y_center + width, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )

        elif x_center > length - radius and y_center < radius:
            # location 4_4
            new_inclusion = self._first_new_inclusion(x_center, y_center, radius, 4)
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center - length, y_center, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array(
                        [x_center - length, y_center + width, radius, 4]
                    ).reshape((1, 4)),
                )
            )
            new_inclusion = np.vstack(
                (
                    new_inclusion,
                    np.array([x_center, y_center + width, radius, 4]).reshape(
                        (1, 4)
                    ),
                )
            )

        else:
            raise Exception(
                "The location of the original point was wrong!!! \n"
            )

        return new_inclusion

    @staticmethod
    def overlap_check(
        new_inclusion: np.ndarray,
        inclusion_pos: np.ndarray,
        dist_factor: float,
        min_dis_threshold: float,
        inclusion_index: int = 0,
        stage: str = "step_one",
        distance_tol: float = 0.025,
    ) -> int:
        """overlap check between new inclusion and the original ones

        Parameters
        ----------
        new_inclusion : np.ndarray
            new inclusion location
        inclusion_pos : np.ndarray
            original inclusions
        dist_factor : float
            distance factor which used to control the minimum distance
            between to inclusions
        inclusion_index : int, optional
            inclusion index , by default 0
        stage : str, optional
            stage of the algorithm, by default "step_one"

        Returns
        -------
        int
            a flag number (1: overlap, 0: non-overlap)

        """

        inclusion_pos = inclusion_pos.copy()

        if stage == "step_one":
            if min_dis_threshold is None:
                min_dis_threshold = dist_factor * (
                    new_inclusion[0, 2] + inclusion_pos[:, 2]
                ).reshape((-1, 1))
            points_dis_temp = distance_matrix(
                inclusion_pos[:, 0:2], new_inclusion[:, 0:2]
            )
            points_dis = np.min(points_dis_temp, 1, keepdims=True)
            min_dis = points_dis - min_dis_threshold

        elif stage == "step_two":
            # calculate the minimum distance threshold
            if min_dis_threshold is None:
                min_dis_threshold = dist_factor * (
                    new_inclusion[0, 2] + inclusion_pos[:, 2]
                ).reshape((-1, 1))
            points_dis_temp = distance_matrix(
                inclusion_pos[:, 0:2], new_inclusion[:, 0:2]
            )
            points_dis_temp[
                inclusion_index: inclusion_index + int(inclusion_pos[inclusion_index, 3]), :
            ] = math.inf
            points_dis = np.min(points_dis_temp, 1, keepdims=True)
            min_dis = points_dis - min_dis_threshold

        else:
            raise ValueError(" Not defined stage \n")

        if min_dis.min() < distance_tol:
            status = 1
        else:
            status = 0

        return status

    @staticmethod
    def min_dis_index(
        temp_inclusion: np.ndarray,
        inclusion_pos: np.ndarray,
        inclusion_min_dis_vector: np.ndarray,
        ii: int,
        cycle: int,
    ) -> list[np.ndarray, int, float]:
        """This function is used to identify the index of closest inclusion
        of every inclusion, which is very import for the first heuristic
        stirring to get more space placing the new disks.

        Parameters
        ----------
        temp_inclusion : np.ndarray
            the inclusion been processed
        inclusion_pos : np.ndarray
            inclusion position
        inclusion_min_dis_vector : np.ndarray
            the first column is the index of the closed point,
            the second column contains the minimum distance between
            those two points
        ii : int
            the index of the being processed point
        cycle : int
            the cycle of the algorithm. At different cycles, the
            criteria of identifying the closed point are  different.

        Returns
        -------
        inclusion_min_dis_vector: np.ndarray
            The updated minimum distance array
        min_index: int
            The index of the minimum distance point
        min_dist : float
            The minimum distance to the minimum distance point

        """
        inclusion_pos = inclusion_pos.copy()

        # pre-process the data : find out the same row data and delete it
        temp_inclusion = temp_inclusion.reshape((1, 2))
        points_dis = distance_matrix(inclusion_pos[:, 0:2], temp_inclusion)
        points_dis[points_dis == 0] = math.inf
        if cycle == 0:
            min_dis = points_dis.min()
            min_index = np.where(points_dis == min_dis)[0]
            inclusion_min_dis_vector[ii, cycle, 0] = min_index
            inclusion_min_dis_vector[ii, cycle, 1] = min_dis
        elif cycle == 1:
            index_pre = int(inclusion_min_dis_vector[ii, cycle - 1, 0])
            if index_pre < points_dis.shape[0]:
                points_dis[index_pre, :] = math.inf
            # identify the minimum index
            min_dis = points_dis.min()
            min_index = np.where(points_dis == min_dis)[0]
            inclusion_min_dis_vector[ii, cycle, 0] = min_index
            inclusion_min_dis_vector[ii, cycle, 1] = min_dis
        else:

            index_pre = int(inclusion_min_dis_vector[ii, cycle - 1, 0])
            index_pre_pre = int(inclusion_min_dis_vector[ii, cycle - 2, 0])
            if (
                index_pre < points_dis.shape[0]
                and index_pre_pre < points_dis.shape[0]
            ):
                points_dis[index_pre, :] = math.inf
                points_dis[index_pre_pre, :] = math.inf
            # identify the minimum index
            min_dis = points_dis.min()
            min_index = np.where(points_dis == min_dis)[0]
            inclusion_min_dis_vector[ii, cycle, 0] = min_index
            inclusion_min_dis_vector[ii, cycle, 1] = min_dis

        return inclusion_min_dis_vector, min_index, min_dis

    @staticmethod
    def gen_heuristic_inclusions(
        ref_point: np.ndarray,
        inclusion_temp: np.ndarray,
        dist_factor: float,
        rng: Any,
    ) -> np.ndarray:
        """Move inclusion to its reference point

        Parameters
        ----------
        ref_point : np.ndarray
            Reference point that the inclusion should move to
        inclusion_temp : np.ndarray
            The considering inclusion
        dist_factor : float
            the minimum distance factor between two inclusions
        rng: Any
            random generator

        Returns
        -------
        np.ndarray
            The updated location of the considering inclusion
        """

        inclusion_temp = inclusion_temp.reshape((1, 3))
        ref_point = ref_point.reshape((1, 3))
        # generate the random factor for inclusion stirring
        delta = rng.uniform(0, 1, 1)
        dist_min = dist_factor * (inclusion_temp[0, 2] + ref_point[0, 2])
        inclusion_loc = inclusion_temp[0, 0:2].reshape((1, 2)).copy()
        ref_loc = ref_point[0, 0:2].reshape((1, 2)).copy()
        # maximum length of movement
        k = 1 - dist_min / distance_matrix(ref_loc, inclusion_loc)
        inclusion_temp[0, 0:2] = inclusion_loc + delta * k * (ref_loc - inclusion_loc)

        return inclusion_temp

    @staticmethod
    def _first_new_inclusion(x: float,
                         y: float,
                         r: float,
                         portion: int) -> np.ndarray:
        """generate the first new inclusion

        Parameters
        ----------
        x : float
            x center
        y : float
            y center
        r : float
            radius
        portion : int
            portion of the inclusion(1, 2, 4)

        Returns
        -------
        np.ndarray
            new inclusion
        """
        return np.array([[x, y, r, portion]])


    @staticmethod
    def inclusion_volume(radius: float) -> float:
        """calculate the inclusion volume of the current inclusion

        Parameters
        ----------
        radius : float
            radius

        Returns
        -------
        vol:float
            volume of current inclusion(disk)
        """
        return np.pi * radius**2

    @staticmethod
    def generate_random_inclusions(
        len_start: float,
        len_end: float,
        wid_start: float,
        wid_end: float,
        radius_mu: float,
        radius_std: float,
        rng,
    ) -> np.ndarray:
        """generate random inclusions with different radiis

        Parameters
        ----------
        len_start : float
            the start location of length
        len_end : float
            the end location of length
        wid_start : float
            the start location of width
        wid_end : float
            the end location of width
        radius_mu : float
            mean of radius
        radius_std : float
            standard deviation of radius
        rng: any
            random seed or generator

        Returns
        -------
        np.ndarray
            location information of generated inclusion
        """

        x = rng.uniform(len_start, len_end, 1)
        y = rng.uniform(wid_start, wid_end, 1)
        r = rng.normal(radius_mu, radius_std, 1)
        # the radius is too small for mesh
        while r <= 0.02*(len_end - len_start - 2*radius_mu):
            r = rng.normal(radius_mu, radius_std, 1)
        inclusion = np.array([x, y, r])
        return inclusion
