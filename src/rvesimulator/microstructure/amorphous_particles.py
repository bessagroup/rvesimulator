import numpy as np

from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation
from typing import List, Tuple
import time

from .base import MicrostructureGenerator


class AmorphousParticles(MicrostructureGenerator):

    def __init__(self,
                 length: int,
                 width: int,
                 mesh_partition_length: int,
                 mesh_partition_width: int,
                 vol_req: List[float],
                 sigma: float,
                 ):
        """ Amorphous particles microstructure generator"""
        self.length = length
        self.width = width
        self.mesh_partition_length = mesh_partition_length
        self.mesh_partition_width = mesh_partition_width
        # define the shape of the image
        self.shape = (mesh_partition_length, mesh_partition_width)
        self.vol_req = vol_req
        self.vol_total = np.sum(vol_req)
        self.sigma = sigma
        # check if the volume fractions
        if len(vol_req) == 1:
            print("Particles with vf: ", vol_req[0])
        elif len(vol_req) == 2:
            print("Particles with vf: ", vol_req[0])
            print("Compatibilizer volume fraction: ", vol_req[1])
        else:
            print("Recheck the volume fractions")

    def generate_microstructure(self, seed=None):

        # fix seed
        np.random.seed(seed)
        # for the first step, generate vf with the desired total vf
        # the matrix has element number of zeros
        # the fiber phase is represented by 1 if no compatibilizer
        # the compatibilizer phase is 1 if there is compatibilizer
        # then the fiber phase is represented by 2
        start_time = time.time()
        # matrix is ones and fibers is zeros
        self.pattern = self._amorphous_shape()
        end_time = time.time()
        self.time_usage = end_time - start_time

    def crate_rgmsh(self):

        # save the microstructure in crate format (matrix should be 1-indexed)
        # fiber should be 2
        crate_discretization = self.pattern.copy()
        crate_discretization[self.pattern == 0] = 2
        self.rgmsh = crate_discretization

    def process_compatibilizer_phase(self,
                                     pattern: np.ndarray) -> Tuple:
        # TODO : check the scheme with Ivan
        # handle the compatibilizer phase by dilating the particles
        structure = np.ones((3, 3))
        dilated_pattern = binary_dilation(pattern,
                                          structure=structure,
                                          border_value=1)
        # add the dilated pattern to the binary pattern
        final_pattern = pattern + dilated_pattern
        # pbc for the compatibilizer phase
        index_interface_horz_edge = np.where(final_pattern[0] == 1)
        # Periodic top-bottom edges
        final_pattern[0][index_interface_horz_edge] = 0
        final_pattern[-1] = final_pattern[0]
        # Periodic left-right edges for the compatibilizer phase
        for i in range(self.shape[1]):
            if final_pattern[i][0] == 1:
                final_pattern[i][0] = 0
            final_pattern[i][-1] = final_pattern[i][0]

        comp_elements = np.where(pattern.flatten() == 1)[0]+1
        pp_elements = np.where(pattern.flatten() == 0)[0]+1
        pe_elements = np.where(pattern.flatten() == 2)[0]+1
        # Verify the compatibilizer volume fraction
        n_compatibilizer_elements = len(comp_elements)
        total_elements = self.mesh_partition_width * \
            self.mesh_partition_length
        if self.vol_req[1] < n_compatibilizer_elements/total_elements:
            idx_selected = np.random.choice(len(comp_elements),
                                            int(self.vol_req[1]
                                                * total_elements),
                                            replace=True)
            idx2remove = np.setdiff1d(
                np.arange(len(comp_elements)), idx_selected)
            #
            pp_elements = list(pp_elements) + \
                list(comp_elements[idx2remove])
            comp_elements = comp_elements[idx_selected]
            return list(pp_elements), list(pe_elements), list(comp_elements)

    def plot_microstructure(self,
                            save_figure=False,
                            fig_name="RVE.png"):

        plt.imshow(self.pattern,
                   cmap='summer_r',
                   origin="lower",
                   interpolation="None",)
        plt.axis('off')
        plt.title(f"vf: {np.sum(self.vol_req)}")
        if save_figure:
            plt.savefig(fig_name, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def to_abaqus_format(self, file_name="micro_structure_info.json"):

        pass

    def gaussian_filter(self,
                        shape: Tuple,
                        sigma: float) -> np.ndarray:
        """ gaussian filter function 

        Parameters
        ----------
        shape : Tuple
            The shape of the filter
        sigma : float
            gaussian filter parameter

        Returns
        -------
        np.ndarray
            The gaussian filter
        """
        x = np.linspace(-shape[1]//2, shape[1]//2, shape[1])
        y = np.linspace(-shape[0]//2, shape[0]//2, shape[0])
        X, Y = np.meshgrid(x, y)
        d = np.sqrt((X / shape[1])**2 + (Y / shape[0])**2)
        return np.exp(-(d**2) / (2 * sigma**2))

    def _amorphous_shape(self):
        """ Generate the amorphous shape with the desired volume fraction
        """
        # generate a random field that fills random values between 0 and 1
        random_pattern = np.random.rand(*self.shape)
        # Fourier Transform into frequency domain
        f_transform = np.fft.fftshift(np.fft.fft2(random_pattern))
        # Gaussian filter in frequency domain
        filter_gaussian = self.gaussian_filter(self.shape, self.sigma)
        f_filtered = f_transform * filter_gaussian
        # Inverse Fourier Transform
        pattern = np.fft.ifft2(np.fft.ifftshift(f_filtered)).real
        # Binary pattern
        threshold = brentq(self.compute_volume_fraction,
                           np.min(pattern),
                           np.max(pattern),
                           args=(pattern, self.vol_total))
        binary_pattern = (pattern < threshold).astype(int)

        # ================== Periodic Boundary Conditions ================== #
        # Periodic top-bottom edges
        binary_pattern[-1] = binary_pattern[0]
        # Periodic top-bottom edges
        binary_pattern[:, -1] = binary_pattern[:, 0]

        return binary_pattern

    @staticmethod
    def compute_volume_fraction(
            threshold: float,
            pattern: np.ndarray,
            desired_volume_fraction: float) -> float:
        """ threshold the pattern to get the desired volume fraction

        Parameters
        ----------
        threshold : float
            threshold value to binarize the pattern
        pattern : np.ndarray
            the random pattern to be threshold
        desired_volume_fraction : float
            the desired volume fraction

        Returns
        -------
        float
            the difference between the desired volume fraction and the
            computed volume fraction
        """
        binary_image = pattern > threshold
        return np.mean(binary_image) - desired_volume_fraction

    @property
    def vol_frac(self):
        return np.sum(self.vol_req)
