#                                                                       Modules
# =============================================================================

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class Simulator:
    """Base class for a FEM simulator"""

    def pre_process(self) -> None:
        """Function that handles the pre-processing"""

        pass

    def execute(self) -> None:
        """Function that calls the FEM simulator the pre-processing"""

        raise NotImplementedError("should be implemented in subclass")

    def post_process(self) -> None:
        """Function that handles the post-processing"""

        raise NotImplementedError("should be implemented in subclass")

    def run(self) -> None:
        """run the simulation"""
        self.pre_process()
        self.execute()
        self.post_process()
