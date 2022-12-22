# third-party
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from SALib.sample import sobol_sequence
from scipy.stats import qmc


class Sampler:
    """design of experiments"""

    def __init__(self) -> None:
        """initilaization"""
        self.seed = None
        self.design_space = None
        self.num_samples = None
        self.num_dim = None
        self.samples = None
        self.responses = None
        self.out_names = None

    def set_seed(self, seed: int = 123456) -> None:
        """This function is used to get the design space info
            and the number of dimension of this problem

            design_space :
        Parameters
        ----------
        seed : int, optional
            a dict that contain the design space, by default 123456
        """

        self.seed = seed

    def set_design_space(self, design_space: dict = None) -> None:
        """This function is used to get the design space info
            and the number of dimension of this problem

        Parameters
        ----------
        design_space : dict, optional
            design space, by default None
        """

        self.design_space = design_space
        self.num_dim = len(design_space)

    def get_samples(self, num_samples: int = None) -> None:
        """ "
            This function is used to get the samples
        Arg:
            num_samples : a dict that contain the design space

        Note:  the function should be completed at the subsclass
        """

        raise NotImplementedError(
            "This function should be implimented in the sub-class \n"
        )

    def create_pandas_frame(self, out_names: list) -> None:
        """
        this function is used to create pandas framework for the doe
        the output will be added at the end of the pandas dataframe but\
        without giving names

        Returns
        -------
        None
        """
        # load the number of outputs and corresponding names
        self.num_outs = len(out_names)

        # transfer the variables to a pandas dataframe
        self.samples = pd.DataFrame(
            self.samples, columns=list(self.design_space.keys())
        )
        responses = np.empty(
            (
                self.num_samples,
                self.num_outs,
            )
        )
        responses[:] = np.nan
        self.responses = pd.DataFrame(responses, columns=out_names)
        self.responses[out_names] = self.responses[out_names].astype(object)

    def save_doe(self, name: str = "doe") -> None:
        """
        This function is used to save the DoE to Json files
        Returns
        -------

        """
        self.samples.to_json(name + ".json", index=True)

    def plot_samples(
        self, fig_name: str = None, save_fig: bool = False
    ) -> None:
        """visualize the two-dimensional sampling


        Parameters
        ----------
        fig_name : str, optional
            figure name , by default None
        save_fig : bool, optional
            sigure figure, by default False
        """
        if self.num_dim == 1:
            with plt.style.context(
                ["ieee", "science", "high-contrast", "grid"]
            ):
                fig, ax = plt.subplots()
                ax.plot(
                    self.samples.values[:, 0],
                    np.zeros((self.num_samples, 1)),
                    ".",
                    label="Samples",
                )
                ax.legend()
                ax.set(xlabel=r"$x_1$")
                ax.set(ylabel=r"$y$")
                ax.autoscale(tight=True)
                if save_fig is True:
                    fig.savefig(fig_name, dpi=300)
                plt.show()

        elif self.num_dim == 2:
            with plt.style.context(
                ["ieee", "science", "high-contrast", "grid"]
            ):
                fig, ax = plt.subplots()
                ax.plot(
                    self.samples.values[:, 0],
                    self.samples.values[:, 1],
                    "*",
                    label="Samples",
                )
                ax.legend()
                ax.set(xlabel=r"$x_1$")
                ax.set(ylabel=r"$x_2$")
                ax.autoscale(tight=True)
                if save_fig is True:
                    fig.savefig(fig_name, dpi=300)
                plt.show()
        else:
            raise Exception("expect 1or2 dimesion problems\n")

    @property
    def samples_(self) -> pd.DataFrame:
        return self.samples

    @property
    def responses_(self) -> pd.DataFrame:
        return self.responses

    @property
    def data(self) -> dict:
        return {"samples": self.samples, "responses": self.responses}


class LatinHyperCube(Sampler):
    def sampling(
        self, num_samples: int, design_space: dict, seed: int, out_names: dict
    ) -> pd.DataFrame:
        self.set_seed(seed=seed)
        self.set_design_space(design_space=design_space)
        self.get_samples(num_samples=num_samples)
        self.create_pandas_frame(out_names=out_names)

        return self.samples

    def get_samples(self, num_samples: int = None) -> np.ndarray:
        self.num_samples = num_samples
        sampler = qmc.LatinHypercube(d=self.num_dim)
        sample = sampler.random(self.num_samples)
        for i, bounds in enumerate(self.design_space.values()):
            sample[:, i] = sample[:, i] * (bounds[1] - bounds[0]) + bounds[0]

        self.samples = sample


class SobolSequence(Sampler):
    def sampling(
        self, num_samples: int, design_space: dict, seed: int, out_names: dict
    ) -> pd.DataFrame:
        self.set_seed(seed=seed)
        self.set_design_space(design_space=design_space)
        self.get_samples(num_samples=num_samples)
        self.create_pandas_frame(out_names=out_names)

        return self.samples

    def get_samples(self, num_samples: int = None) -> np.ndarray:
        self.num_samples = num_samples
        sample = sobol_sequence.sample(self.num_samples, self.num_dim)
        for i, bounds in enumerate(self.design_space.values()):
            sample[:, i] = sample[:, i] * (bounds[1] - bounds[0]) + bounds[0]
        self.samples = sample


class RandomSampler(Sampler):
    def sampling(
        self, num_samples: int, design_space: dict, seed: int, out_names: dict
    ) -> pd.DataFrame:
        self.set_seed(seed=seed)
        self.set_design_space(design_space=design_space)
        self.get_samples(num_samples=num_samples)
        self.create_pandas_frame(out_names=out_names)

        return self.samples

    def get_samples(self, num_samples: int = None) -> np.ndarray:
        self.num_samples = num_samples
        np.random.seed(self.seed)
        sample = np.random.random((self.num_samples, self.num_dim))
        for i, bounds in enumerate(self.design_space.values()):
            sample[:, i] = sample[:, i] * (bounds[1] - bounds[0]) + bounds[0]

        self.samples = sample


class FixNumberSampler(Sampler):
    def sampling(
        self, num_samples: int, design_space: dict, out_names: dict
    ) -> pd.DataFrame:
        self.set_design_space(design_space=design_space)
        self.get_samples(num_samples=num_samples)
        self.create_pandas_frame(out_names=out_names)
        return self.samples

    def get_samples(self, num_samples: int = None) -> np.ndarray:
        self.num_samples = num_samples
        fixedvalue = list(self.design_space.values())
        sample = np.repeat(fixedvalue[0], num_samples)
        self.samples = sample
