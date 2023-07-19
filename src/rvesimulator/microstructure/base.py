#                                                                       Modules
# =============================================================================

# Third party

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


class MicrosctuctureGenerator:
    "base class of mirostructure generator"

    def generate_microstructure(self, seed: any = None) -> float:
        """generating micro-structure

        Parameters
        ----------
        seed : any, optional
            seed generator or number , by default None
        """

        raise NotImplementedError("The function should be implemented \
                                  in sub-class \n")

    def plot_microstructure(
        self, save_figure: bool = False, fig_name: str = "RVE.png"
    ) -> None:
        """plot figure for RVE

        Parameters
        ----------
        save_figure : bool, optional
            save figure or not , by default False
        fig_name : str, optional
            figure name, by default "RVE.png"

        Raises
        ------
        NotImplementedError
            error report
        """

        raise NotImplementedError("The function should be implemented \
                                  in sub-class \n")

    def to_abaqus_format(self,
                         file_name: str = "micro_structure_info.json") -> None:
        """convert microstructure to abaqus format

        Parameters
        ----------
        file_name : str, optional
            file name, by default "micro_structure_info.json"

        Raises
        ------
        NotImplementedError
            error report
        """

        raise NotImplementedError("The function should be implemented \
                                  in sub-class \n")

    def to_crate_format(self, file_name: str = "microstructure.rgmsh") -> None:
        """convert microstructure to crate format

        Parameters
        ----------
        file_name : str, optional
            _description_, by default "microstructure.rgmsh"
        """

        np.save(file_name, self.rgmsh.T)

    def rgmsh_plot(
        self,
        save_fig: bool = False,
        fig_name: str = "rgmsh.png",
        **kwarg,
    ) -> None:
        """contour or voxel plot of the discrete microstructure

        Parameters
        ----------
        save_fig : bool, optional
            save figure, by default False
        fig_name : str, optional
            figure name, by default "rgmsh.png"
        """
        if self.fiber_positions.shape[1] == 4:
            fig = plt.figure(**kwarg)
            ax = fig.add_subplot(111)
            ax.imshow(
                self.rgmsh.T,
                origin="lower",
                interpolation="None",
                cmap="summer_r",
            )
            ax.set_yticks([])
            ax.set_xticks([])
            if save_fig is True:
                plt.savefig(fig_name, dpi=300)
                plt.close()
            plt.show()
        elif self.fiber_positions.shape[1] == 5:
            fig = plt.figure(**kwarg)
            ax = fig.add_subplot(projection="3d")
            ax.voxels(self.rgmsh.T, facecolors="#EE3377")
            ax.set_yticks([])
            ax.set_xticks([])
            if save_fig is True:
                plt.savefig(fig_name, dpi=300)
                plt.close()
            plt.show()

    @staticmethod
    def circle_plot(
        fibers: np.ndarray,
        length: float,
        width: float,
        vol_frac: float,
        save_figure: bool = False,
        fig_name: str = "microstrcuture.png",
        **kwargs,
    ) -> None:
        """2d rve plot for circle inclusion

        Parameters
        ----------
        fibers : np.ndarray
            fiber locations/positions
        length : float
            length of the rve
        width : float
            with of the rve
        vol_frac : float
            reached volume fraction
        save_figure : bool, optional
            save figure , by default False
        fig_name : str, optional
            figure name , by default "microstrcuture.png"
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
        plt.xlim((0.0, length))
        plt.ylim((0.0, width))
        plt.title(f"$V_f$ = {vol_frac*100:.2f}")
        axes.set_yticks([])
        axes.set_xticks([])
        if save_figure is True:
            plt.savefig(fig_name, dpi=300)
            plt.close()
        plt.show()

    @staticmethod
    def sphere_coordinate(location_information: np.ndarray) -> tuple:
        """generate coordinate of sphere

        Parameters
        ----------
        location_information : np.ndarray
           a numpy array contains center and radius info of
           a sphere [x, y, z, r]

        Returns
        -------ss
        tuple
            coordinate of x, y, z
        """

        loc_info = np.reshape(location_information, (1, -1))
        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 20)
        x = loc_info[0, 3] * np.outer(np.cos(u), np.sin(v)) + loc_info[0, 0]
        y = loc_info[0, 3] * np.outer(np.sin(u), np.sin(v)) + loc_info[0, 1]
        z = loc_info[0, 3] * np.outer(np.ones(np.size(u)), np.cos(v)) \
            + loc_info[0, 2]

        return x, y, z

    def sphere_plot(
        self,
        fibers: np.ndarray,
        length: float,
        width: float,
        height: float,
        vol_frac: float,
        save_figure: bool = False,
        fig_name: str = "cubic.png",
    ) -> None:
        """plot 3d rve with sphere inclusion

        Parameters
        ----------
        location_information : np.ndarray
            location information
        length : float
            length of cubic
        width : float
            width of cubic
        height : float
            height of cubic
        vol_frac : float
            volume fraction
        save_figure : bool, optional
            a flag to indicate if we want save the figure or not,
            by default False
        fig_name : str, optional
            fig name, by default "cubic_rve.png"
        """
        # TODO refine this function with better color scheme
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        # Plot the surface
        for ii in range(fibers.shape[0]):
            x, y, z = \
                self.sphere_coordinate(location_information=fibers[ii, :])
            if fibers[ii, 4] == 1:
                ax.plot_surface(x, y, z, color="lightseagreen")
            elif fibers[ii, 4] == 2:
                ax.plot_surface(x, y, z, color="orange")
            elif fibers[ii, 4] == 4:
                ax.plot_surface(x, y, z, color="blue")
            elif fibers[ii, 4] == 8:
                ax.plot_surface(x, y, z, color="firebrick")
            else:
                raise ValueError("Spltting number is wrong!\n")
        axes = [int(length), int(width), int(height)]
        data = np.ones(axes, dtype=np.bool_)
        alpha = 0.6
        colors = np.empty(axes + [4], dtype=np.float32)
        colors[:] = [1, 1, 1, alpha]
        ax.voxels(data, facecolors=colors)
        ax.set_xlim3d(0, length)
        ax.set_ylim3d(0, width)
        ax.set_zlim3d(0, height)
        plt.title(f"$V_f$ = {vol_frac*100:.2f}")

        # Set an equal aspect ratio
        ax.set_aspect("auto")
        if not save_figure:
            plt.show()
        else:
            plt.savefig(fig_name, dpi=300, bbox_inches="tight")
            plt.close()
