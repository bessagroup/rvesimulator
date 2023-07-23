from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class HardeningLaw:
    """hardening law"""

    hardening_law_table: np.ndarray = None

    def linear(self, **kwargs) -> list:
        """linear hardening law

        Returns
        -------
        hardening_law_table : list
            hardening law table
        """
        # assert the needed constants
        assert "a" in kwargs.keys(), "provide 'a' value "
        assert "yield_stress" in kwargs.keys(), "provide 'yield_stress' value"
        # get the arguements
        yield_stress = kwargs["yield_stress"]
        a = kwargs["a"]
        # dfine the table for hardening law
        hardening_law_table = np.zeros((101, 2))
        hardening_law_table[:, 1] = np.linspace(0, 1, 101)

        # generate the hardening law
        hardening_law_table[:, 0] = (
            yield_stress + a * hardening_law_table[:, 1]
        )
        hardening_law_table[-1, 1] = 10.0
        hardening_law_table[-1, 0] = (
            yield_stress + a * hardening_law_table[-1, 1]
        )

        hardening_law_table = hardening_law_table.T
        self.hardening_law_table = hardening_law_table
        return hardening_law_table.tolist()

    def swift(self, **kwargs) -> list:
        """swift hardening law

        Returns
        -------
        hardening_law_table : list
            hardening law table
        """
        # assert the needed constants
        assert "a" in kwargs.keys(), "provide 'a' value "
        assert "b" in kwargs.keys(), "provide 'b' value "
        assert "yield_stress" in kwargs.keys(), "provide 'yield_stress' value"
        # get the arguements
        yield_stress = kwargs["yield_stress"]
        a = kwargs["a"]
        b = kwargs["b"]
        # dfine the table for hardening law
        hardening_law_table = np.zeros((101, 2))
        hardening_law_table[:, 1] = np.linspace(0, 1, 101)
        # generate the hardening law
        hardening_law_table[:, 0] = (
            yield_stress + a * (hardening_law_table[:, 1]) ** b
        )
        hardening_law_table[-1, 1] = 10.0
        hardening_law_table[-1, 0] = (
            yield_stress + a * (hardening_law_table[-1, 1]) ** b
        )
        hardening_law_table = hardening_law_table.T
        self.hardening_law_table = hardening_law_table
        return hardening_law_table.tolist()

    def ramberg(self, **kwargs) -> list:
        """ramberg hardening law

        Returns
        -------
        hardening_law_table : list
            hardening law table
        """
        # assert the needed constants
        assert "a" in kwargs.keys(), "provide 'a' value "
        assert "b" in kwargs.keys(), "provide 'b' value "
        assert "yield_stress" in kwargs.keys(), "provide 'yield_stress' value"
        # get the arguments
        yield_stress = kwargs["yield_stress"]
        a = kwargs["a"]
        b = kwargs["b"]
        # define the table for hardening law
        hardening_law_table = np.zeros((101, 2))
        hardening_law_table[:, 1] = np.linspace(0, 1, 101)

        # generate the hardening law
        hardening_law_table[:, 0] = yield_stress * (
            1 + a * (hardening_law_table[:, 1])
        ) ** (1 / b)
        hardening_law_table[-1, 1] = 10.0
        hardening_law_table[-1, 0] = yield_stress * (
            1 + a * (hardening_law_table[-1, 1])
        ) ** (1 / b)

        hardening_law_table = hardening_law_table.T
        self.hardening_law_table = hardening_law_table

        return hardening_law_table.tolist()

    def hardening_law_plot(self, **kwargs) -> None:
        "plot the hardening law"
        fig, ax = plt.subplots(**kwargs)
        ax.plot(
            self.hardening_law_table[1, :],
            self.hardening_law_table[0, :],
            label="hardening_law",
        )
        plt.legend()
        plt.xlabel(r"$\varepsilon$")
        plt.ylabel(r"$\sigma$")
        plt.xlim([0, 1])
        plt.grid("-.")
        plt.show()
