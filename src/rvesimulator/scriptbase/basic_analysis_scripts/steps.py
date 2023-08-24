# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *


class FieldOutputs:
    """define field outputs
    """

    def _define_field_outputs(self):
        # create Final-outputs
        self.model.fieldOutputRequests["F-Output-1"].setValues(
            variables=(
                "S",
                "E",
                "LE",
                "ENER",
                "ELEN",
                "ELEDEN",
                "EVOL",
                "IVOL",
            ),
            timeInterval=self.time_interval,
        )
        self.model.FieldOutputRequest(
            name="F-Output-2",
            createStepName="Step-1",
            variables=("U", "RF"),
            timeInterval=self.time_interval,
        )


class HistoryOutputs:
    """define history outputs
    """
    def _define_history_outputs(self):
        self.model.historyOutputRequests["H-Output-1"].setValues(
            variables=(
                "ALLAE",
                "ALLCD",
                "ALLIE",
                "ALLKE",
                "ALLPD",
                "ALLSE",
                "ALLWK",
                "ETOTAL",
            ),
            timeInterval=self.time_interval,
        )


class SmallDeformationSteps(FieldOutputs, HistoryOutputs):
    """define small deformation steps

    Parameters
    ----------
    FieldOutputs : class
        field outputs
    HistoryOutputs : class
        history outputs
    """
    def create_step(
        self,
        initialInc=0.1,
        maxInc=1,
        maxNumInc=1000000,
        minInc=1e-20,
    ):

        self.model.StaticStep(name="Step-1", previous="Initial")
        self.model.StaticStep(
            initialInc=initialInc,
            maxInc=maxInc,
            maxNumInc=maxNumInc,
            minInc=minInc,
            name="Step-1",
            previous="Initial",
            timePeriod=self.time_period,
        )
        # define the outputs
        self._define_history_outputs()
        self._define_field_outputs()


class LargeDeformationSteps(FieldOutputs, HistoryOutputs):
    """large deformation steps

    Parameters
    ----------
    FieldOutputs : class
        filed outputs
    HistoryOutputs : class
        historical outputs
    """
    def create_step(
        self,
        initialInc=0.1,
        maxInc=1,
        maxNumInc=1000000,
        minInc=1e-20,
    ):

        self.model.StaticStep(name="Step-1", previous="Initial")
        self.model.StaticStep(
            initialInc=initialInc,
            maxInc=maxInc,
            maxNumInc=maxNumInc,
            minInc=minInc,
            name="Step-1",
            previous="Initial",
            timePeriod=self.time_period,
            nlgeom=ON,
        )

        # define the outputs
        self._define_history_outputs()
        self._define_field_outputs()


