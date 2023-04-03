import os

import matplotlib.pyplot as plt
import numpy as np

from rvesimulator.microstructure.circle_particles import CircleParticles
from rvesimulator.microstructure.shpere_particles import SphereParticles


class CircleMircoStructure:
    def __call__(
        self,
        size: float,
        radius_mu: float,
        radius_std: float,
        vol_req: float,
        seed: any,
    ) -> None:

        self.microstructure_generator = CircleParticles(
            length=size,
            width=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
        )
        self.microstructure_generator.generate_microstructure(seed=seed)
        microstrcutrue_info = self.microstructure_generator.to_abaqus_format(
            file_name="micro_structure.json"
        )
        self.microstructure_generator.plot_microstructure(save_figure=True)

        return microstrcutrue_info, self.microstructure_generator.vol_frac


class SphereMicroStructrue:
    def __call__(
        self,
        size: float,
        radius_mu: float,
        radius_std: float,
        vol_req: float,
        seed: any = None,
    ) -> None:
        self.microstructure_generator = SphereParticles(
            length=size,
            width=size,
            height=size,
            radius_mu=radius_mu,
            radius_std=radius_std,
            vol_req=vol_req,
        )
        self.microstructure_generator.generate_microstructure(seed=seed)
        microstrcutrue_info = self.microstructure_generator.to_abaqus_format(
            file_name="micro_structure.json"
        )
        return microstrcutrue_info, self.microstructure_generator.vol_frac


class CircleSVEMicroStructure:
    def __init__(self, size: float = 1.0, radius: float = 0.2) -> None:
        self.size = size
        self.radius = radius

    def location_information(self, task: str = "task1") -> np.ndarray:

        size = self.size
        radius = self.radius
        # compatibility check
        if radius > 0.17 * size:
            raise Exception(
                "The radius should smaller than 1/8 of the size \n"
            )
        # generate fiber location
        if task == "task1":
            location_information = [
                [size / 2, size / 4, radius, 1],
                [size / 2, 3 * size / 4, radius, 1],
            ]
        elif task == "task2":
            location_information = [
                [size / 4, size / 4, radius, 1],
                [size / 4, 3 * size / 4, radius, 1],
                [3 * size / 4, size / 4, radius, 1],
                [3 * size / 4, 3 * size / 4, radius, 1],
            ]
        elif task == "task3":
            location_information = [
                [size / 2, size / 4, 2 * radius, 1],
                [size / 2, 3 * size / 4, 2 * radius, 1],
            ]
        elif task == "task4":
            location_information = [
                [size / 4, size / 4, 2 * radius, 1],
                [size / 4, 3 * size / 4, 2 * radius, 1],
                [3 * size / 4, size / 4, 2 * radius, 1],
                [3 * size / 4, 3 * size / 4, 2 * radius, 1],
                [size / 2, size / 2, radius, 1],
            ]
        return location_information


def cricle_plot(
    fibers: np.ndarray,
    length: float,
    width: float,
    **kwargs,
) -> None:

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
    axes.set_yticks([])
    axes.set_xticks([])
    plt.show()
