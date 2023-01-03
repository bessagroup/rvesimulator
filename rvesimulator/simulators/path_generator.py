import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


class PathGenerator:
    def __init__(
        self, num_control_points: int = 6, num_increment: int = 100
    ) -> None:
        """initialization of path generation

        Parameters
        ----------
        num_control_points : int, optional
            number of control point, by default 6
        """
        self.num_control_point = num_control_points
        self.num_increment = num_increment

        # generate control points
        self.x_control = np.linspace(
            0, num_increment, num_control_points, endpoint=True
        )
        self.y_control = np.zeros((1, 3))
        self.y_control = np.vstack(
            (
                self.y_control,
                np.random.uniform(-1, 1, (num_control_points - 1, 3)),
            )
        )

        # define the output raviables
        self.path = None

    def linear_interpolate(self) -> None:
        """linear interpolate"""

        fit1 = interp1d(self.x_control, self.y_control[:, 0], kind="linear")
        fit2 = interp1d(self.x_control, self.y_control[:, 1], kind="linear")
        fit3 = interp1d(self.x_control, self.y_control[:, 2], kind="linear")

        # test points
        self.x_test = np.linspace(
            0, self.num_increment, self.num_increment + 1, endpoint=True
        )
        self.y1 = fit1(self.x_test)
        self.y2 = fit2(self.x_test)
        self.y3 = fit3(self.x_test)

        self.path = np.array([self.y1, self.y2, self.y3]).tolist()

        return self.path

    def quadratic_interpolate(self) -> None:
        """quadratic interpolate"""
        fit1 = interp1d(self.x_control, self.y_control[:, 0], kind="quadratic")
        fit2 = interp1d(self.x_control, self.y_control[:, 1], kind="quadratic")
        fit3 = interp1d(self.x_control, self.y_control[:, 2], kind="quadratic")

        # test points
        self.x_test = np.linspace(
            0, self.num_increment, self.num_increment + 1, endpoint=True
        )

        self.y1 = fit1(self.x_test)
        self.y2 = fit2(self.x_test)
        self.y3 = fit3(self.x_test)

        self.path = np.array([self.y1, self.y2, self.y3]).tolist()

        return self.path

    def plot_path(
        self, fig_name: str = "loads_path.png", save_fig: bool = False
    ) -> None:
        """plot the loads path"""

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
        ax.plot(self.x_test, self.y1, label="path xx")
        ax.plot(self.x_test, self.y2, label="path yy")
        ax.plot(self.x_test, self.y3, label="path xy")
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
