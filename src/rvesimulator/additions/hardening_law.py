"""
Material hardening laws for Abaqus default material models
"""

#                                                                       Modules
# =============================================================================

# Standard
from abc import ABC
from typing import List

# Third party
import matplotlib.pyplot as plt
import numpy as np

#                                                        Authorship and Credits
# =============================================================================
__author__ = 'Jiaxiang Yi (J.Yi@tudelft.nl)'
__credits__ = ['Jiaxiang Yi']
__status__ = 'Stable'
# =============================================================================
#
# =============================================================================


#                                                                       Classes
# =============================================================================

class HardeningLaw(ABC):
    def __init__(self) -> None:
        """Abstract class for the hardening law. The hardening law is defined
        by the hardening law table. The hardening law table is a 2D array
        with the first row being the strain and the second row being the
        stress.
        """

    def calculate_hardening_table(self) -> List[float]:
        ...

    def plot_hardening_law(self, **kwargs):
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
        plt.grid("-.")
        plt.show()

#                                                                Hardening laws
# =============================================================================


class LinearHardeningLaw(HardeningLaw):
    def __init__(self,
                 a: float = 0.2,
                 yield_stress: float = 0.5) -> None:
        self.a = a
        self.yield_stress = yield_stress

    def calculate_hardening_table(self) -> List[float]:
        # get the arguments
        yield_stress = self.yield_stress
        a = self.a
        # define the table for hardening law
        hardening_law_table = np.zeros((201, 2))
        hardening_law_table[:, 1] = np.linspace(0, 2, 201)

        # generate the hardening law
        hardening_law_table[:, 0] = (
            yield_stress + a * hardening_law_table[:, 1]
        )
        hardening_law_table[-1, 1] = 10.0
        hardening_law_table[-1, 0] = (
            yield_stress + a * hardening_law_table[-1, 1]
        )

        self.hardening_law_table = hardening_law_table.T
        return self.hardening_law_table.tolist()


class SwiftHardeningLaw(HardeningLaw):
    def __init__(self,
                 a: float = 0.2,
                 b: float = 0.2,
                 yield_stress: float = 0.5) -> None:
        self.a = a
        self.b = b
        self.yield_stress = yield_stress

    def calculate_hardening_table(self) -> List[float]:
        # get the arguments
        yield_stress = self.yield_stress
        a = self.a
        b = self.b
        eps = 1e-5
        # define the table for hardening law
        hardening_law_table = np.zeros((201, 2))
        hardening_law_table[:, 1] = np.linspace(0, 2, 201)
        # generate the hardening law
        hardening_law_table[:, 0] = (
            yield_stress + a * (hardening_law_table[:, 1]) ** b
        )
        hardening_law_table[-1, 1] = 10.0
        hardening_law_table[-1, 0] = (
            yield_stress + a * (hardening_law_table[-1, 1] + eps) ** b
        )
        self.hardening_law_table = hardening_law_table.T
        return self.hardening_law_table.tolist()


class RambergHardeningLaw(HardeningLaw):

    def __init__(self,
                 a: float = 0.2,
                 b: float = 0.2,
                 yield_stress: float = 0.5) -> None:
        self.a = a
        self.b = b
        self.yield_stress = yield_stress

    def calculate_hardening_table(self) -> List[float]:
        # get the arguments
        yield_stress = self.yield_stress
        a = self.a
        b = self.b
        # define the table for hardening law
        hardening_law_table = np.zeros((201, 2))
        hardening_law_table[:, 1] = np.linspace(0, 2, 201)

        # generate the hardening law
        hardening_law_table[:, 0] = yield_stress * (
            1 + a * (hardening_law_table[:, 1])
        ) ** (1 / b)
        hardening_law_table[-1, 1] = 10.0
        hardening_law_table[-1, 0] = yield_stress * (
            1 + a * (hardening_law_table[-1, 1])
        ) ** (1 / b)

        self.hardening_law_table = hardening_law_table.T
        return self.hardening_law_table.tolist()
