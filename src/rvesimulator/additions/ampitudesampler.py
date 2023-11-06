#                                                                       Modules
# =============================================================================
# Standard
from typing import Any, Tuple

# Third-party
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

#                                                          Authorship & Credits
# =============================================================================
__author__ = 'Jiaxiang Yi'
__credits__ = ['Jiaxiang Yi']
__status__ = 'Stable'
# =============================================================================

class AmplitudeGenerator:
    """amplitude generator"""

    def __init__(self, num_dim: int = 3) -> None:
        """Initialization"""
        self.num_dim = num_dim

        # assert dimension
        assert num_dim == 3 or num_dim == 6, "dimension should be 3 or 6"

    def get_amplitude(
        self,
        num_amplitude: int,
        num_control: int,
        num_steps: int,
        arg_name: str = "amplitude",
        interpolation_method: str = "quadratic",
        seed: Any = None,
    ) -> pd.DataFrame:
        """get amplitude curves

        Parameters
        ----------
        num_amplitude : int
            number of amplitude, usually 3 or 6
        num_control : int
            control points number
        num_steps : int
            num steps of ABAQUS simulation
        arg_name : str, optional
            name of amplitude, by default "amplitude"
        interpolation_method : str, optional
            interpolation method, by default "quadratic"
        seed : any, optional
            seed , by default None

        Returns
        -------
        pd.DataFrame
            _description_
        """
        self.arg_name = arg_name
        self.num_steps = num_steps
        # create df Dataframe for amplitude
        amplitude_temp = np.empty([num_amplitude, 3])
        amplitude_temp[:] = np.nan
        self.amplitude = pd.DataFrame(
            amplitude_temp, columns=["x_control", "y_control", arg_name]
        ).astype(object)
        # create pandas frame for output
        amplitude_to_abaqus = np.empty([num_amplitude, 1])
        self.amplitude_to_abaqus = pd.DataFrame(
            amplitude_to_abaqus, columns=[arg_name]
        ).astype(object)

        # a loop to generate loads path
        for ii in range(num_amplitude):
            # generate control points
            x_control, y_control = self.generate_control_points(
                seed=seed + ii,
                num_control=num_control,
                num_steps=num_steps,
                num_dim=self.num_dim,
            )
            self.amplitude.at[ii, "x_control"] = x_control
            self.amplitude.at[ii, "y_control"] = y_control
            # save for plot figure
            self.amplitude.at[ii, arg_name] = self.interpolation(
                x_control=x_control,
                y_control=y_control,
                num_steps=num_steps,
                interpolation_method=interpolation_method,
            )

            # save for abaqus simulation
            self.amplitude_to_abaqus.at[ii, arg_name] = self.amplitude.at[
                ii, arg_name
            ].T.tolist()

        return self.amplitude_to_abaqus

    @staticmethod
    def generate_control_points(
        seed: int, num_control: int, num_steps: int, num_dim: int = 3
    ) -> Tuple[Any, Any]:
        """generate control points

        Parameters
        ----------
        seed : int
            seed
        num_control : int
            control points number
        num_steps : int
            number of steps for abaqus simulation
        num_dim : int, optional
            number of dimension, by default 3

        Returns
        -------
        Tuple[Any, Any]
            x_location of control points, y_location of control points
        """
        num_control = int(num_control)
        num_increment = int(num_steps)
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
        num_steps: int,
        interpolation_method: str = "quadratic",
    ) -> np.ndarray:
        """do interpolation

        Parameters
        ----------
        x_control : np.ndarray
            x location of control points
        y_control : np.ndarray
            y location of control points
        num_steps : int
            number of steps of abaqus simulation
        interpolation_method : str, optional
            interpolation method, by default "quadratic"

        Returns
        -------
        np.ndarray
            amplitude
        """
        num_dim = y_control.shape[1]
        num_steps = int(num_steps)
        path = np.zeros((num_steps + 1, num_dim))
        for ii in range(num_dim):
            fit = interp1d(
                x_control, y_control[:, ii], kind=interpolation_method
            )
            # x_increment points
            num_increment = np.linspace(
                0, num_steps, num_steps + 1, endpoint=True
            )
            # get the path
            path[:, ii] = fit(num_increment)

            # clear the fit
            del fit

        return path

    def plot_amplitude(
        self,
        fig_name: str = "amplitude_path.png",
        save_fig: bool = False,
        internal: bool = True,
        **kwargs,
    ) -> None:
        """plot amplitude

        Parameters
        ----------
        fig_name : str, optional
            figure name, by default "amplitude_path.png"
        save_fig : bool, optional
            save figure, by default False
        internal : bool, optional
            approach this function from internal or external, by default True
        """
        if internal is True:
            assert (
                "iteration" in kwargs.keys()
            ), " iteration should be provided"
            iteration = kwargs["iteration"]
            x_control = self.amplitude.at[iteration, "x_control"]
            y_control = self.amplitude.at[iteration, "y_control"]
            amplitude = self.amplitude.at[iteration, self.arg_name]
            num_steps = self.num_steps
        else:
            assert (
                "x_control" in kwargs.keys()
            ), " x_control should be provided"
            assert (
                "y_control" in kwargs.keys()
            ), " y_control should be provided"
            assert (
                "num_steps" in kwargs.keys()
            ), " num_steps should be provided"
            assert "amplitude" in kwargs.keys(), "  path should be provided"
            x_control = kwargs["x_control"]
            y_control = kwargs["y_control"]
            num_steps = kwargs["num_steps"]
            amplitude = kwargs["amplitude"]
        num_steps = int(num_steps)
        x_increment = np.linspace(0, num_steps, num_steps + 1, endpoint=True)
        # two dimension or three dimension problem (from the view of FEM)
        if y_control.shape[1] == 3:
            fig, ax = plt.subplots()
            ax.plot(x_control, y_control[:, 0], "o", label="xx points")
            ax.plot(x_control, y_control[:, 1], "o", label="yy points")
            ax.plot(x_control, y_control[:, 2], "o", label="xy points")
            ax.plot(x_increment, amplitude[:, 0], label="xx path")
            ax.plot(x_increment, amplitude[:, 1], label="yy path")
            ax.plot(x_increment, amplitude[:, 2], label="xy path")
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
            fig, ax = plt.subplots()
            ax.plot(x_control, y_control[:, 0], "o", label="xx points")
            ax.plot(x_control, y_control[:, 1], "o", label="yy points")
            ax.plot(x_control, y_control[:, 2], "o", label="zz points")
            ax.plot(x_control, y_control[:, 3], "o", label="xy points")
            ax.plot(x_control, y_control[:, 4], "o", label="xz points")
            ax.plot(x_control, y_control[:, 5], "o", label="yz points")
            ax.plot(x_increment, amplitude[:, 0], label="path xx")
            ax.plot(x_increment, amplitude[:, 1], label="path yy")
            ax.plot(x_increment, amplitude[:, 2], label="path zz")
            ax.plot(x_increment, amplitude[:, 3], label="path xy")
            ax.plot(x_increment, amplitude[:, 4], label="path xz")
            ax.plot(x_increment, amplitude[:, 5], label="path yz")
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
