#                                                                       Modules
# =============================================================================
# Third party
import numpy as np
import pandas as pd
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


class StrainPathSampler:
    def __init__(self, num_dim: int = 3, seed: int = 123) -> None:
        """Initialization of strain path sampler"""

        self.data = None
        self.path = None
        self.seed = seed
        self.num_dim = num_dim

        # assert the dimension
        assert num_dim == 3 or num_dim == 6, "dimension should be 3 or 6"

    def get_strain_path(
        self,
        data: dict = None,
        arg_name: str = "loads_path",
        interploation_method: str = "quadratic",
    ) -> dict:
        """get the strain path

        Parameters
        ----------
        data : dict, optional
            samples with number of control point and number of increments,
            by default None

        Returns
        -------
        data:dict
            updated data dict with a new colume named strain path
        """

        # get the data dicts
        self.data = data.copy()
        self.num_points = len(data["samples"].axes[0])
        self.argument = arg_name
        # create an empty dataframe that contain the keys of loads_path
        path_temp = np.empty([self.num_points, 3])
        path_temp[:] = np.nan
        self.path = pd.DataFrame(
            path_temp, columns=["x_control", "y_control", arg_name]
        )
        self.path[["x_control", "y_control", arg_name]] = self.path[
            ["x_control", "y_control", arg_name]
        ].astype(object)
        # create df dataframe for output
        loads_path_temp = np.empty([self.num_points, 1])
        loads_path_temp[:] = np.nan
        loads_path = pd.DataFrame(loads_path_temp, columns=[arg_name])
        loads_path[arg_name] = loads_path[arg_name].astype(object)
        # a loop to generate loads path
        for ii in range(self.num_points):
            # generate control points
            x_control, y_control = self.generate_control_points(
                seed=self.seed + ii,
                num_control=data["samples"].at[ii, "num_control"],
                num_increment=data["samples"].at[ii, "num_increment"],
                num_dim=self.num_dim,
            )
            # save for plot figure
            self.path.at[ii, "x_control"] = x_control
            self.path.at[ii, "y_control"] = y_control
            self.path.at[ii, arg_name] = self.interpolation(
                x_control=x_control,
                y_control=y_control,
                num_increment=data["samples"].at[ii, "num_increment"],
                interploation_method=interploation_method,
            )
            # save for abaqus simulation
            loads_path.at[ii, arg_name] = self.path.at[ii, arg_name].tolist()
        self.data["samples"] = pd.concat(
            [data["samples"], loads_path], axis=1, join="inner"
        )

        return self.data

    @staticmethod
    def generate_control_points(
        seed: int, num_control: int, num_increment: int, num_dim: int = 3
    ) -> tuple[
        np.ndarray[any, np.dtype[np.floating]],
        np.ndarray[any, np.dtype[np.floating]],
    ]:
        """generate control points

        Parameters
        ----------
        seed : int
            seed for repeating
        num_control : int
            num of control points
        num_increment : int
            number of increment

        Returns
        -------
        x_control, y_control: tuple[ np.ndarray[Any, np.dtype[np.floating]],
        np.ndarray[Any, np.dtype[np.floating]], ]
            control points for x and y axises
        """
        num_control = np.int64(num_control)
        num_increment = np.int64(num_increment)
        np.random.seed(seed)
        x_control = np.linspace(0, num_increment, num_control, endpoint=True)
        y_control = np.zeros((1, num_dim))
        y_control = np.vstack(
            (
                y_control,
                np.random.uniform(-1, 1, (num_control - 1, num_dim)),
            )
        )
        return x_control, y_control

    @staticmethod
    def interpolation(
        x_control: np.ndarray,
        y_control: np.ndarray,
        num_increment: int,
        interploation_method: str = "quadratic",
    ) -> np.ndarray:
        """_summary_

        Parameters
        ----------
        x_control : np.ndarray
            _description_
        y_control : np.ndarray
            _description_
        num_increment : int
            _description_
        interploation_method : str, optional
            _description_, by default "quadratic"

        Returns
        -------
        _type_
            _description_
        """
        # identify the dimension
        num_dim = y_control.shape[1]
        num_increment = np.int64(num_increment)
        path = np.zeros((num_increment + 1, num_dim))
        for ii in range(num_dim):
            fit = interp1d(
                x_control, y_control[:, ii], kind=interploation_method
            )
            # x_increment points
            x_increment = np.linspace(
                0, num_increment, num_increment + 1, endpoint=True
            )
            # get the path
            path[:, ii] = fit(x_increment)

            # clear the fit
            del fit

        return path

    def plot_path(
        self,
        fig_name: str = "loads_path.png",
        save_fig: bool = False,
        internal: bool = True,
        **kwarg,
    ) -> None:

        """plot path

        Parameters
        ----------
        fig_name : str, optional
            name of the figure, by default "loads_path.png"
        save_fig : bool, optional
            save figure or not, by default False
        """
        # enough arguments should be probided
        if internal is True:
            assert "iteration" in kwarg.keys(), " iteration should be provided"
            iteration = kwarg["iteration"]
            x_control = self.path.at[iteration, "x_control"]
            y_control = self.path.at[iteration, "y_control"]
            num_increment = self.data["samples"].at[iteration, "num_increment"]
            path = self.path.at[iteration, self.argument]
        else:
            assert "x_control" in kwarg.keys(), " x_control should be provided"
            assert "y_control" in kwarg.keys(), " y_control should be provided"
            assert (
                "num_increment" in kwarg.keys()
            ), " num_increment should be provided"
            assert "path" in kwarg.keys(), "  path should be provided"
            x_control = kwarg["x_control"]
            y_control = kwarg["y_control"]
            num_increment = kwarg["num_increment"]
            path = kwarg["path"]
        num_increment = np.int64(num_increment)
        x_increment = np.linspace(
            0, num_increment, num_increment + 1, endpoint=True
        )
        # two dimension or three dimension problem (from the view of FEM)
        if y_control.shape[1] == 3:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.plot(x_control, y_control[:, 0], "o", label="xx points")
            ax.plot(x_control, y_control[:, 1], "o", label="yy points")
            ax.plot(x_control, y_control[:, 2], "o", label="xy points")
            ax.plot(x_increment, path[:, 0], label="xx path")
            ax.plot(x_increment, path[:, 1], label="yy path")
            ax.plot(x_increment, path[:, 2], label="xy path")
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
        elif y_control.shape[1] == 6:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.plot(x_control, y_control[:, 0], "o", label="xx points")
            ax.plot(x_control, y_control[:, 1], "o", label="yy points")
            ax.plot(x_control, y_control[:, 2], "o", label="zz points")
            ax.plot(x_control, y_control[:, 3], "o", label="xy points")
            ax.plot(x_control, y_control[:, 4], "o", label="xz points")
            ax.plot(x_control, y_control[:, 5], "o", label="yz points")
            ax.plot(x_increment, path[:, 0], label="path xx")
            ax.plot(x_increment, path[:, 1], label="path yy")
            ax.plot(x_increment, path[:, 2], label="path zz")
            ax.plot(x_increment, path[:, 3], label="path xy")
            ax.plot(x_increment, path[:, 4], label="path xz")
            ax.plot(x_increment, path[:, 5], label="path yz")
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
