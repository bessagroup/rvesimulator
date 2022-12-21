# system packages
import json
import math
import time
from math import fmod

# third-party packages
import numpy as np
from scipy.spatial import distance_matrix

# import local functions
from rvesimulator.microstructures.microstrcuture_2d import (
    MicrosctucturaGenerator2D,
)
from rvesimulator.microstructures.microstructure_plots import DrawRVE2D


class CircleInclusion(MicrosctucturaGenerator2D, DrawRVE2D):
    """
    This class is used to generate the 2D RVE with disks particles

    """

    def __init__(
        self,
        length: float,
        width: float,
        radius: float,
        vol_req: float,
        num_guess_max: int = 50000,
        num_fiber_max: int = 750,
        num_cycle_max: int = 15,
        dist_min_factor: float = 2.07,
        second_heuristic: bool = True,
    ) -> None:

        """Initialization

        Parameters
        ----------
        length : float
            length of RVE
        width : float
            width of RVE
        radius : float
            radius of fibers
        vol_req : float
            required volume fraction
        num_guess_max : int, optional
            maximum guess for fibers, by default 50000
        num_fiber_max : int, optional
            maximum fiber inside RVE, by default 750
        num_cycle_max : int, optional
            iteration cycles, by default 15
        dist_min_factor : float, optional
            distance factor, by default 2.07
        second_heuristic: bool, optional
            second stage heuristic, by default True
        """

        # geometry information of the 2D RVE with homogeneous circles
        self.length = length
        self.width = width
        self.radius = radius
        self.vol_req = vol_req

        # Initialization of the algorithm
        self.dist_min = dist_min_factor * self.radius
        # Number of guesses in the random generation of fibre position
        self.num_guess_max = num_guess_max
        self.num_fibers_max = num_fiber_max
        self.num_cycles_max = num_cycle_max

        self.second_heuristic = second_heuristic

    def __parameter_initialization(self) -> None:
        """Initialize the parameters"""
        self.num_change = 3
        self.num_cycle = 0
        self.vol_frac = 0
        self.vol_per_fiber = np.pi * self.radius**2
        # Total volume of image
        self.vol_total = self.length * self.width
        # initial size of square in Second Heuristic
        self.square_size = 3 * self.radius
        self.square_inc = (8.5 - 10 * self.vol_req) * self.radius

        # initial coordinate for position of fibre
        self.len_start = -1 * self.radius
        self.len_end = self.length + self.radius
        self.wid_start = -1 * self.radius
        self.wid_end = self.width + self.radius

        # some import variables
        self.fiber_positions = None

    def generate_rve(self) -> float:
        """core function for generating RVE microstructures via Melro's
        algorithm

        Returns
        -------
        float
            Actual volume fracture
        """

        start_time = time.time()
        self.__generate_rve()
        # check the periodical compatability
        while np.count_nonzero(self.fiber_positions[:, 2] == 2) % 2 != 0:
            print(
                f"the half fiber numbers \
                {np.count_nonzero(self.fiber_positions[:, 2] == 2)} \n"
            )

            self.__generate_rve()
        end_time = time.time()
        print(
            f"Time of generate the 2D RVE with volume fraction\
            = {self.vol_frac*100:.2f} $%$ is {(end_time - start_time):.2f} s"
        )

        return self.vol_frac

    def save_results(
        self, file_name: str = "micro_structure_info.json"
    ) -> None:
        """save the needed information of the 2D \
            microstructure into "Json" file.
        Later this Json file will be imported to ABAQUS
        """

        geometry_info = {
            "location_information": self.fiber_positions.tolist(),
            "radius": self.radius,
            "len_start": self.len_start,
            "wid_start": self.wid_start,
            "len_end": self.len_end,
            "wid_end": self.wid_end,
        }
        with open(file_name, "w") as fp:
            json.dump(geometry_info, fp)

    def plot_rve(
        self, save_figure: bool = False, fig_name: str = "RVE.png"
    ) -> None:
        """
        This function is used to visualize the 2D RVE
        """
        self.cricle_inclusion_plot(
            circle_position=self.fiber_positions,
            radius=self.radius,
            len_start=self.len_start,
            len_end=self.len_end,
            wid_start=self.wid_start,
            wid_end=self.wid_end,
            vol_frac=self.vol_frac,
            save_figure=save_figure,
            fig_name=fig_name,
        )

    def __generate_rve(self) -> None:
        """core procedure of generating RVE"""
        self.__parameter_initialization()
        self._procedure_initialization()
        self._core_iteration()

    def _procedure_initialization(self) -> None:
        """This function is used to generate the first disk and
        assign the initial values of the algorithm
        """

        # initialization(generate the first fiber randomly)
        self.fiber_min_dis_vector = np.zeros(
            (self.num_fibers_max, self.num_cycles_max + 1, 2)
        )
        self.num_fibers = 1
        # generate the location of the first fiber
        fiber_temp = self.generate_random_fibers(
            len_start=self.radius,
            len_end=self.length - self.radius,
            wid_start=self.radius,
            wid_end=self.width - self.radius,
        )
        # update the volume fraction information
        self.vol = self.vol_per_fiber / self.vol_total
        self.fiber_positions = np.zeros((1, 3))
        self.fiber_positions[0, 0] = fiber_temp[0, 0]
        self.fiber_positions[0, 1] = fiber_temp[1, 0]
        self.fiber_positions[0, 2] = 1

    def _core_iteration(self) -> None:
        """core iteration part of the micro-structure generation method"""

        # the main loop of the algorithm
        while (
            self.vol_frac < self.vol_req
            and self.num_cycle < self.num_cycles_max
        ):
            ############################################################
            ########        generate the inclusions randomly   #########
            ############################################################
            self.num_trial = 1
            while (
                self.num_trial < self.num_guess_max
                and self.vol_frac < self.vol_req
                and self.num_fibers < self.num_fibers_max
            ):
                # update the info of number trial
                self.num_trial = self.num_trial + 1
                fiber_temp = self.generate_random_fibers(
                    len_start=self.len_start,
                    len_end=self.len_end,
                    wid_start=self.wid_start,
                    wid_end=self.wid_end,
                )
                # check the location of the fiber and
                new_fiber = self.new_positions(
                    x_center=fiber_temp[0, 0],
                    y_center=fiber_temp[1, 0],
                    radius=self.radius,
                    length=self.length,
                    width=self.width,
                )
                # check the overlap of new fiber
                overlap_status = self.overlap_check(
                    new_fiber=new_fiber,
                    fiber_pos=self.fiber_positions,
                    dist_min=self.dist_min,
                )
                if overlap_status == 0:
                    self.fiber_positions = np.vstack(
                        (self.fiber_positions, new_fiber)
                    )
                    self.vol_frac = (
                        self.vol_frac + self.vol_per_fiber / self.vol_total
                    )
                    self.num_fibers = self.num_fibers + new_fiber.shape[0]
                del new_fiber

            ############################################################
            ##### First Heuristic -Move Fibers to gain more areas ######
            ############################################################
            ii = 0
            if self.fiber_positions.shape[0] < self.num_fibers_max:
                # for every point, stirring is needed!!!
                # print('Begin first heuristic stirring \n')
                while ii < self.fiber_positions.shape[0]:
                    (
                        self.fiber_min_dis_vector,
                        min_index,
                        min_dis,
                    ) = self.min_dis_index(
                        self.fiber_positions[ii, 0:2],
                        self.fiber_positions,
                        self.fiber_min_dis_vector,
                        ii,
                        self.num_cycle,
                    )
                    new_fiber_temp = self.generate_first_heuristic_fibers(
                        ref_point=self.fiber_positions[min_index, 0:2],
                        fiber_temp=self.fiber_positions[ii, 0:2],
                        dist_min=self.dist_min,
                    )
                    new_fiber = self.new_positions(
                        x_center=new_fiber_temp[0, 0],
                        y_center=new_fiber_temp[0, 1],
                        radius=self.radius,
                        length=self.length,
                        width=self.width,
                    )
                    overlap_status = self.overlap_check(
                        new_fiber=new_fiber,
                        fiber_pos=self.fiber_positions,
                        dist_min=self.dist_min,
                        stage="step_two",
                        fiber_index=ii,
                    )
                    # check: if the new fibers(cause it maybe more than
                    # 1 fiber centers) will overlap with the
                    # remaining ones or not
                    if overlap_status == 0:
                        ii = self._update_fiber_position(
                            new_fiber=new_fiber, iter=ii
                        )
                    else:
                        ii = ii + int(self.fiber_positions[ii, 2])
                    del new_fiber, new_fiber_temp

            ############################################################
            ###Second Heuristic-Move Fibers around to gain more areas ##
            ############################################################
            if self.second_heuristic is True:
                # TODO check the correctness of this stage
                if (self.fiber_positions.shape[0] < self.num_fibers_max) and (
                    (fmod(self.num_cycle, 2) == 0)
                    or (self.square_size > (self.length / 2 - self.radius))
                ):
                    jj = 0
                    while jj < self.fiber_positions.shape[0]:
                        # identify the region of the point
                        fiber_temp = self.fiber_positions[jj, :]
                        # find the new location of this fiber
                        fiber_new_pos = self._second_heuristic_stirring(
                            fiber_temp, jj
                        )
                        if fiber_new_pos is None:
                            jj = jj + int(fiber_temp[2])
                        else:
                            new_fiber = self.new_positions(
                                x_center=fiber_new_pos[0, 0],
                                y_center=fiber_new_pos[0, 1],
                                radius=self.radius,
                                length=self.length,
                                width=self.width,
                            )
                            # it will generate a vector (split = 1, 2, 4)
                            # verify the overlap
                            overlap_status = self.overlap_check(
                                new_fiber=new_fiber,
                                fiber_pos=self.fiber_positions,
                                dist_min=self.dist_min,
                                stage="step_three",
                                fiber_index=jj,
                            )
                            # remaining ones or not
                            if overlap_status == 0:
                                jj = self._update_fiber_position(
                                    new_fiber=new_fiber, iter=jj
                                )
                            else:
                                jj = jj + int(fiber_temp[2])
                            del new_fiber
                        del fiber_new_pos, fiber_temp
                    # update the square size
                    if self.square_size <= (self.length / 2 - self.radius):
                        self.square_size = self.square_size + self.square_inc
            # end of one cycle
            self.num_cycle = self.num_cycle + 1

    def _update_fiber_position(self, new_fiber: np.ndarray, iter: int) -> int:
        """update the fiber position

        Parameters
        ----------
        new_fiber : np.ndarray
            the generated new fiber
        iter : int
            determine the nex fiber should be analysis

        Returns
        -------
        iter : int
            he updated number of index
        """
        # check the location compatibility
        if new_fiber[0, 2] != new_fiber.shape[0]:
            raise Exception("fiber number is wrong \n")

        if self.fiber_positions[iter, 2] == 1:
            # the original fiber is a full circle
            if new_fiber[0, 2] == 1:
                # keep one point after stirring
                self.fiber_positions[iter, :] = new_fiber
                iter = iter + 1
            elif new_fiber[0, 2] == 2:
                # split into two hal after stirring
                self.fiber_positions = np.delete(
                    self.fiber_positions, (iter), axis=0
                )
                # delete the original ii point first
                self.fiber_positions = np.insert(
                    self.fiber_positions, (iter, iter + 1), [new_fiber], axis=0
                )
                # add two point at iith location
                iter = iter + 2
            elif new_fiber[0, 2] == 4:
                # split into foul one-fouth circles after stirring
                self.fiber_positions = np.delete(
                    self.fiber_positions, (iter), axis=0
                )
                # delete the original ii point first
                self.fiber_positions = np.insert(
                    self.fiber_positions,
                    (iter, iter + 1, iter + 2, iter + 3),
                    [new_fiber],
                    axis=0,
                )
                # add four point at iith location
                iter = iter + 4
            else:
                raise KeyError("The new splits of fiber is wrong!!! \n")

        elif self.fiber_positions[iter, 2] == 2:
            if new_fiber[0, 2] == 2:
                # keep two half points after stirring
                self.fiber_positions[iter : iter + 2, :] = new_fiber
                iter = iter + 2
            elif new_fiber[0, 2] == 1:
                # reduce to one circle after stirring
                self.fiber_positions = np.delete(
                    self.fiber_positions, (iter, iter + 1), axis=0
                )
                # delete the ii+1 half fiber
                self.fiber_positions = np.insert(
                    self.fiber_positions, (iter), [new_fiber], axis=0
                )
                # add four point at iith location
                iter = iter + 1
            elif new_fiber[0, 2] == 4:
                # reduce to one circle after stirring
                self.fiber_positions = np.delete(
                    self.fiber_positions, (iter, iter + 1), axis=0
                )
                # delete the ii+1 half fiber
                self.fiber_positions = np.insert(
                    self.fiber_positions,
                    (iter, iter + 1, iter + 2, iter + 3),
                    [new_fiber],
                    axis=0,
                )
                # add four point at iith location
                iter = iter + 4
            else:
                print("The new splits of fiber is wrong!!! \n")

        elif self.fiber_positions[iter, 2] == 4:
            if new_fiber[0, 2] == 4:
                # reduce to one circle after stirring
                self.fiber_positions[iter : iter + 4, :] = new_fiber
                iter = iter + 4
            elif new_fiber[0, 2] == 2:
                # reduce to one circle after stirring
                self.fiber_positions = np.delete(
                    self.fiber_positions,
                    (iter, iter + 1, iter + 2, iter + 3),
                    axis=0,
                )
                # delete the original ii point first
                self.fiber_positions = np.insert(
                    self.fiber_positions, (iter, iter + 1), [new_fiber], axis=0
                )
                # add two point at iith location
                iter = iter + 2
            elif new_fiber[0, 2] == 1:
                # reduce to one circle after stirring
                self.fiber_positions = np.delete(
                    self.fiber_positions,
                    (iter, iter + 1, iter + 2, iter + 3),
                    axis=0,
                )
                # delete the original ii point first
                self.fiber_positions = np.insert(
                    self.fiber_positions, (iter), [new_fiber], axis=0
                )
                # add two point at iith location
                iter = iter + 1
            else:
                raise KeyError("The new splits of fiber is wrong!!! \n")

        else:
            raise KeyError("Can not match the overlap condition  !!! \n")
            iter = iter + int(self.fiber_positions[iter, 2])

        return iter

    def _region_identifier(self, fiber_position: np.ndarray) -> np.ndarray:
        """this function aims to predict which location the disk been
        process belongs to
        ####################################
        #     #          2_1       #       #
        # 2_0  #                    #2_2   #
        ####################################
        #     #                    #       #
        # 1_0 #        1_1         #   1_2 #
        #     #                    #       #
        #     #                    #       #
        ####################################
        # 0_0 #         0_1        # 0_2   #
        #     #                    #       #
        ####################################


        Parameters
        ----------
        fiber_position : np.ndarray
            the location of the been processed fiber

        Returns
        -------
        index_matrix:np.ndarray
             a location information of the been processed disk
        """

        # note first index in the box is index y and the second
        fiber_position = fiber_position.copy()
        fiber_position = fiber_position.reshape((1, 3))
        index_x = np.digitize(
            fiber_position[0, 0],
            bins=np.array([self.square_size, self.width - self.square_size]),
        )
        index_y = np.digitize(
            fiber_position[0, 1],
            bins=np.array([self.square_size, self.length - self.square_size]),
        )
        index_matrix = np.array([index_y, index_x])
        index_matrix = index_matrix.reshape(1, 2)

        return index_matrix

    def _second_heuristic_stirring(
        self, fiber_position: np.ndarray, iter: int
    ) -> np.ndarray:
        """core part of second heuristic stirring

        Parameters
        ----------
        fiber_position : np.ndarray
            the disk which is being processed
        iter : int
            the index of the disk been processed

        Returns
        -------
        fiber_new_position: np.ndarray
            the updated location of the disk
        """

        start_angel = np.array(
            [
                [0.0, 0.0, np.pi / 2],
                [-np.pi / 2, None, np.pi / 2],
                [3 * np.pi / 2, np.pi, np.pi],
            ]
        )
        end_angel = np.array(
            [
                [np.pi / 2, np.pi, np.pi],
                [np.pi / 2, None, 3 * np.pi / 2],
                [2 * np.pi, 2 * np.pi, 3 * np.pi / 2],
            ]
        )
        rad = np.array(
            [0.75 * self.radius, 0.5 * self.radius, 0.25 * self.radius]
        )
        fiber_position = fiber_position.copy()
        for ii, rad_temp in enumerate(rad):
            region_index = self._region_identifier(
                fiber_position=fiber_position
            )
            # if the point located in center, then it just stay
            # the original location
            if region_index[0, 1] == region_index[0, 0] == 1:
                fiber_new_position = None
                break
            else:
                stir_theta = np.linspace(
                    start_angel[region_index[0, 0], region_index[0, 1]],
                    end_angel[region_index[0, 0], region_index[0, 1]],
                    90,
                )
                fiber_new_positions = self.generate_second_heuristic_fibers(
                    fiber_position, rad_temp, stir_theta
                )
                # check the overlap info and return the points meet the
                # stirring criteria
                fiber_new_position = self.overlap_check_second_heuristic(
                    fiber_new_positions,
                    self.fiber_positions,
                    self.dist_min,
                    iter,
                )
                if fiber_new_position is None:
                    continue
                else:
                    break

        return fiber_new_position

    @staticmethod
    def generate_random_fibers(
        len_start: float, len_end: float, wid_start: float, wid_end: float
    ) -> np.ndarray:
        """_summary_

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

        Returns
        -------
        np.ndarray
            location information of generated fiber
        """

        x = np.random.uniform(len_start, len_end, 1)
        y = np.random.uniform(wid_start, wid_end, 1)
        fiber = np.array([x, y])
        return fiber

    @staticmethod
    def generate_first_heuristic_fibers(
        ref_point: np.ndarray, fiber_temp: np.ndarray, dist_min: float
    ) -> np.ndarray:
        """Move fiber to its reference point

        Parameters
        ----------
        ref_point : np.ndarray
            Reference point that the fiber should move to
        fiber_temp : np.ndarray
            The considering fiber
        dist_min : float
            the minimum distance requirement between two fibers

        Returns
        -------
        np.ndarray
            The updated location of the considering fiber
        """

        fiber_temp = fiber_temp.reshape((1, 2))
        # generate the random factor for fiber stirring
        delta = np.random.uniform(0, 1, 1)
        # maximum length of movement
        k = 1 - dist_min / distance_matrix(ref_point, fiber_temp)
        fiber_temp = fiber_temp + delta * k * (ref_point - fiber_temp)

        return fiber_temp

    @staticmethod
    def generate_second_heuristic_fibers(
        fiber_temp: np.ndarray, rad: float, stir_theta: float
    ) -> np.ndarray:
        """Generating a sequential points of temp fiber locations

        Parameters
        ----------
        fiber_temp : np.ndarray
            the considering fiber
        rad : float
            the considering angle
        stir_theta : float
            the stirring theta

        Returns
        -------
        np.ndarray
            the updated fiber location
        """

        fiber_temp = fiber_temp.reshape(1, 3)
        fiber_location = np.array(
            [rad * np.cos(stir_theta), rad * np.sin(stir_theta)]
        ).reshape(90, 2)
        fiber_location = fiber_location + fiber_temp[0, 0:2].reshape(-1, 2)

        return fiber_location

    @staticmethod
    def overlap_check_second_heuristic(
        fiber_temp: np.ndarray,
        fibers_location: np.ndarray,
        dist_min: float,
        iter: int,
    ) -> any:
        """Check the if the new fiber will overlap with previous ones in
        second heuristic stage.

        Parameters
        ----------
        fiber_temp : np.ndarray
            location of the fiber (180 rows)
        Fibers_location : np.ndarray
            the existed fiber locations
        dist_min : float
            the allowed minimum distance between fibers
        iter : int
            the index of the fiber being processed

        Returns
        -------
        any
            the new fiber or None
        """

        # check the overlap with the existed points first
        # delete the iter points temperary
        fibers_location_copy = fibers_location.copy()
        fibers_location_copy = np.delete(fibers_location_copy, (iter), axis=0)
        points_dis = distance_matrix(fiber_temp, fibers_location_copy[:, 0:2])
        points_dis[points_dis == 0] = math.inf
        points_dis_min = points_dis.min(axis=1)
        # print(points_dis)
        overlap_matrix = np.digitize(points_dis.min(axis=1), bins=[dist_min])
        if np.sum(overlap_matrix) == 0:
            return None
        else:
            min_index = np.where(overlap_matrix == 1)
            selected_min = points_dis_min[min_index].min()
            min_index_final = np.where(points_dis == selected_min)[0]
            fiber_new = fiber_temp[min_index_final, :]
            return fiber_new

    @staticmethod
    def new_positions(
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
            XC, YC, split  which is the new locations of this fiber
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

        new_fiber = np.zeros((1, 3))
        if (
            radius <= x_center <= length - radius
            and radius <= y_center <= width - radius
        ):
            # locate in center region and split = 1
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 1

        elif length - radius >= x_center >= radius > y_center:
            # location 2_1
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 2
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center, y_center + width, 2]).reshape((1, 3)),
                )
            )

        elif (
            radius <= x_center <= length - radius and y_center > width - radius
        ):
            # location 2_2
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 2
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center, y_center - width, 2]).reshape((1, 3)),
                )
            )

        elif width - radius >= y_center >= radius > x_center:
            # location 2_3
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 2
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center + length, y_center, 2]).reshape((1, 3)),
                )
            )

        elif (
            radius <= y_center <= width - radius and x_center > length - radius
        ):
            # location 2_4
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 2
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center - length, y_center, 2]).reshape((1, 3)),
                )
            )

        elif x_center < radius and y_center > width - radius:
            # location 4_1
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 4
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center + length, y_center, 4]).reshape((1, 3)),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center + length, y_center - width, 4]).reshape(
                        (1, 3)
                    ),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center, y_center - width, 4]).reshape((1, 3)),
                )
            )

        elif x_center > length - radius and y_center > width - radius:
            # location 4_2
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 4
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center - length, y_center, 4]).reshape((1, 3)),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center - length, y_center - width, 4]).reshape(
                        (1, 3)
                    ),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center, y_center - width, 4]).reshape((1, 3)),
                )
            )

        elif x_center < radius and y_center < radius:
            # location 4_3
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 4
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center + length, y_center, 4]).reshape((1, 3)),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center + length, y_center + width, 4]).reshape(
                        (1, 3)
                    ),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center, y_center + width, 4]).reshape((1, 3)),
                )
            )

        elif x_center > length - radius and y_center < radius:
            # location 4_4
            new_fiber[0, 0] = x_center
            new_fiber[0, 1] = y_center
            new_fiber[0, 2] = 4
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center - length, y_center, 4]).reshape((1, 3)),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center - length, y_center + width, 4]).reshape(
                        (1, 3)
                    ),
                )
            )
            new_fiber = np.vstack(
                (
                    new_fiber,
                    np.array([x_center, y_center + width, 4]).reshape((1, 3)),
                )
            )

        else:
            raise Exception(
                "The location of the original point was wrong!!! \n"
            )

        return new_fiber

    @staticmethod
    def overlap_check(
        new_fiber: np.ndarray,
        fiber_pos: np.ndarray,
        dist_min: float,
        fiber_index: int = 0,
        stage: str = "step_one",
    ) -> int:
        """Check if the new positions will overlap with the original one
        or not ?
        Parameters

        Parameters
        ----------
        new_fiber : np.ndarray
            he locations of the new point(in some cases: several
            positions will be generated for one point because of
            the splitting )
        fiber_pos : np.ndarray
            location of the original points
        dist_min : float
            the allowed minimum distance between disks
        index : int, optional
            the index of the checking fibers, which should not be 0
            when excuting stage two and three
        stage : str, optional
            the stage of the algorithm, step_one, step_two,step_three
            , by default "step_one"

        Returns
        -------
        int
            overlap status, 0 is non-overlap and 1 means there is
            overlap
        """
        fiber_pos = fiber_pos.copy()

        if stage == "step_one":
            points_dis = distance_matrix(fiber_pos[:, 0:2], new_fiber[:, 0:2])
            min_dis = points_dis.min()

        elif stage == "step_two" or stage == "step_three":
            points_dis = distance_matrix(fiber_pos[:, 0:2], new_fiber[:, 0:2])
            points_dis[
                fiber_index : fiber_index + int(fiber_pos[fiber_index, 2]), :
            ] = math.inf
            min_dis = points_dis.min()
        else:
            raise ValueError(" Not defined stage \n")

        if min_dis < dist_min:
            status = 1
        else:
            status = 0

        return status

    @staticmethod
    def min_dis_index(
        temp_fiber: np.ndarray,
        fiber_pos: np.ndarray,
        fiber_min_dis_vector: np.ndarray,
        ii: int,
        cycle: int,
    ) -> list[np.ndarray, int, float]:
        """This function is used to identify the index of closest fiber
        of every fiber, which is very import for the first heuristic
        stirring to get more space placing the new disks.

        Parameters
        ----------
        temp_fiber : np.ndarray
            the fiber been processed
        fiber_pos : np.ndarray
            fiber position
        fiber_min_dis_vector : np.ndarray
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
        fiber_min_dis_vector: np.ndarray
            The updated minimum distance array
        two: int
            The index of the minimum distance point
        three : float
            The minimum distance to the minimum distance point

        """
        fiber_pos = fiber_pos.copy()

        # pre-process the data : find out the same row data and delete it
        temp_fiber = temp_fiber.reshape((1, 2))
        points_dis = distance_matrix(fiber_pos[:, 0:2], temp_fiber)
        points_dis[points_dis == 0] = math.inf
        if cycle == 0:
            min_dis = points_dis.min()
            min_index = np.where(points_dis == min_dis)[0]
            fiber_min_dis_vector[ii, cycle, 0] = min_index
            fiber_min_dis_vector[ii, cycle, 1] = min_dis
        elif cycle == 1:
            index_pre = int(fiber_min_dis_vector[ii, cycle - 1, 0])
            if index_pre < points_dis.shape[0]:
                points_dis[index_pre, :] = math.inf
            # identify the minimum index
            min_dis = points_dis.min()
            min_index = np.where(points_dis == min_dis)[0]
            fiber_min_dis_vector[ii, cycle, 0] = min_index
            fiber_min_dis_vector[ii, cycle, 1] = min_dis
        else:

            index_pre = int(fiber_min_dis_vector[ii, cycle - 1, 0])
            index_pre_pre = int(fiber_min_dis_vector[ii, cycle - 2, 0])
            if (
                index_pre < points_dis.shape[0]
                and index_pre_pre < points_dis.shape[0]
            ):
                points_dis[index_pre, :] = math.inf
                points_dis[index_pre_pre, :] = math.inf
            # identify the minimum index
            min_dis = points_dis.min()
            min_index = np.where(points_dis == min_dis)[0]
            fiber_min_dis_vector[ii, cycle, 0] = min_index
            fiber_min_dis_vector[ii, cycle, 1] = min_dis

        return fiber_min_dis_vector, min_index, min_dis
