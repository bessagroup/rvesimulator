#                                                                       Modules
# =============================================================================
# Third party
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class PathGenerator:
    def __init__(
        self, num_control_points: int = 6, num_increment: int = 100
    ) -> None:
        """initialization of path generation

        Parameters
        ----------
        num_control_points : int, optional
            number of control point, by default 6
        num_increment: int, optional
            number of increment of simulation, by default 100
        """
        self.num_control_point = np.int64(num_control_points)
        self.num_increment = np.int64(num_increment)

        # generate control points
        self.x_control = np.linspace(
            0, self.num_increment, self.num_control_point, endpoint=True
        )
        self.y_control = np.zeros((1, 3))
        self.y_control = np.vstack(
            (
                self.y_control,
                np.random.uniform(-1, 1, (self.num_control_point - 1, 3)),
            )
        )

        # define the output raviables
        self.path = None

    def linear_interpolate(self) -> None:
        """linear interpolate"""

        fit1 = interp1d(self.x_control, self.y_control[:, 0], kind="linear")
        fit2 = interp1d(self.x_control, self.y_control[:, 1], kind="linear")
        fit3 = interp1d(self.x_control, self.y_control[:, 2], kind="linear")

        # x_increment points
        self.x_increment = np.linspace(
            0, self.num_increment, self.num_increment + 1, endpoint=True
        )
        self.path1 = fit1(self.x_increment)
        self.path2 = fit2(self.x_increment)
        self.path3 = fit3(self.x_increment)

        self.path = np.array([self.path1, self.path2, self.path3]).tolist()

        return self.path

    def quadratic_interpolate(self) -> None:
        """quadratic interpolate"""
        fit1 = interp1d(self.x_control, self.y_control[:, 0], kind="quadratic")
        fit2 = interp1d(self.x_control, self.y_control[:, 1], kind="quadratic")
        fit3 = interp1d(self.x_control, self.y_control[:, 2], kind="quadratic")

        # x_increment points
        self.x_increment = np.linspace(
            0, self.num_increment, self.num_increment + 1, endpoint=True
        )
        self.path1 = fit1(self.x_increment)
        self.path2 = fit2(self.x_increment)
        self.path3 = fit3(self.x_increment)

        self.path = np.array([self.path1, self.path2, self.path3]).tolist()

        return self.path

    def plot_path(
        self, fig_name: str = "loads_path.png", save_fig: bool = False
    ) -> None:
        """plot path

        Parameters
        ----------
        fig_name : str, optional
            name of the figure, by default "loads_path.png"
        save_fig : bool, optional
            save figure or not, by default False
        """

        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot(
            self.x_control,
            self.y_control[:, 0],
            "o",
            label="path xx control points",
        )
        ax.plot(
            self.x_control,
            self.y_control[:, 1],
            "o",
            label="path yy control points",
        )
        ax.plot(
            self.x_control,
            self.y_control[:, 2],
            "o",
            label="path xy control points",
        )
        ax.plot(self.x_increment, self.path1, label="path xx")
        ax.plot(self.x_increment, self.path2, label="path yy")
        ax.plot(self.x_increment, self.path3, label="path xy")
        plt.legend(
            bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
            loc="lower left",
            mode="expand",
            borderaxespad=0.0,
            ncol=2,
        )
        if save_fig is True:
            fig.savefig(fname=fig_name, dpi=300)
        else:
            plt.show()
