import matplotlib.pyplot as plt
import numpy as np


class DrawRVE2D:
    @staticmethod
    def cricle_inclusion_plot(
        circle_position: np.ndarray,
        radius: float,
        len_start: float,
        len_end: float,
        wid_start: float,
        wid_end: float,
        vol_frac: float,
        save_figure: bool = False,
        fig_name: str = "RVE.png",
    ) -> None:
        """vistualize the 2D RVE

        Parameters
        ----------
        circle_position : np.ndarray
            center position of fibers
        radius : float
            radius of fibers
        len_start : float
            start length of RVE modeling in Abaqus
        len_end : float
            end length of RVE modeling in Abaqus
        wid_start : float
            start width of RVE modeling in Abaqus
        wid_end : float
            end width of RVE modeling in Abaqus
        vol_frac : float
            Actual volume fraction
        save_figure : bool, optional
            save figure , by default False
        """
        # number of inclusions
        num_circles = circle_position.shape[0]
        figure, axes = plt.subplots()
        for ii in range(num_circles):
            if circle_position[ii, 2] == 1:
                cc = plt.Circle(
                    (circle_position[ii, 0], circle_position[ii, 1]),
                    radius,
                    color="lightseagreen",
                )
            elif circle_position[ii, 2] == 2:
                cc = plt.Circle(
                    (circle_position[ii, 0], circle_position[ii, 1]),
                    radius,
                    color="orange",
                )
            elif circle_position[ii, 2] == 4:
                cc = plt.Circle(
                    (circle_position[ii, 0], circle_position[ii, 1]),
                    radius,
                    color="firebrick",
                )
            else:
                raise ValueError("Spltting number is wrong!  \n")
            axes.add_artist(cc)
        axes.set_aspect(1)
        plt.vlines(
            x=len_start + radius,
            ymin=wid_start + radius,
            ymax=wid_end - radius,
            colors="red",
            ls="--",
        )
        plt.vlines(
            x=len_end - radius,
            ymin=wid_start + radius,
            ymax=wid_end - radius,
            colors="red",
            ls="--",
        )
        plt.hlines(
            xmin=len_start + radius,
            xmax=len_end - radius,
            y=wid_end - radius,
            colors="red",
            ls="--",
        )
        plt.hlines(
            xmin=len_start + radius,
            xmax=len_end - radius,
            y=wid_start + radius,
            colors="red",
            ls="--",
        )
        plt.vlines(
            x=len_start, ymin=wid_start, ymax=wid_end, colors="green", ls=":"
        )
        plt.vlines(
            x=len_end, ymin=wid_start, ymax=wid_end, colors="green", ls=":"
        )
        plt.hlines(
            xmin=len_start, xmax=len_end, y=wid_end, colors="green", ls=":"
        )
        plt.hlines(
            xmin=len_start, xmax=len_end, y=wid_start, colors="green", ls=":"
        )
        plt.xlim((len_start - radius, len_end + radius))
        plt.ylim((wid_start - radius, wid_end + radius))
        plt.title(f"$V_f$ = {vol_frac*100:.2f}")
        if not save_figure:
            plt.show()
            plt.close()
        else:
            plt.savefig(fig_name, dpi=300, bbox_inches="tight")
            plt.close()
