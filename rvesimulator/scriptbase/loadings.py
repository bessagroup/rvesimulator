#                                                                       Modules
# =============================================================================
# Third party
import numpy

# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Alpha"
# =============================================================================
#
# =============================================================================


class Loading2D:
    def __init__(self, model, assembly, instance_name):
        """initialization of loadings

        Parameters
        ----------
        model : abaqus model
            abaqus model
        assembly : assembly
            abaqus assembly
        instance_name : instance
            abaqus instance
        """
        self.model = model
        self.assembly = assembly
        self.instance = instance_name

    def create_loads(self, loads):
        # loads should be a list with three elemnts

        if loads[1] == loads[2] == 0.0:
            # x axis tension
            self.strain_11(E11=loads[0])
            self.rigid_cons_11()
        elif loads[0] == loads[1] == 0.0:
            # y axis tension
            self.strain_22(E22=loads[1])
            self.rigid_cons_22()
        elif loads[0] == loads[2] == 0.0:
            # pure shear
            self.strain_12(E12=loads[2])
        else:
            self.complex_loads(loads=loads)

    def strain_11(self, E11):
        self.model.DisplacementBC(
            name="E_11",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=E11,
            u2=UNSET,
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )

    def strain_22(self, E22):
        self.model.DisplacementBC(
            name="E_22",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-T"],
            u1=UNSET,
            u2=E22,
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )

    def strain_12(self, E12):
        self.model.DisplacementBC(
            name="E_12",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=UNSET,
            u2=E12,
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
        self.model.DisplacementBC(
            name="E_21",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-T"],
            u1=E12,
            u2=UNSET,
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )

    def complex_loads(self, loads):
        self.model.DisplacementBC(
            name="E_R_Ref",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=loads[0],
            u2=loads[2],
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
        self.model.DisplacementBC(
            name="E_T_Ref",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-T"],
            u1=loads[2],
            u2=loads[1],
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )

    def rigid_cons_11(self):
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="RigidConsE111",
            region=self.assembly.instances[self.instance].sets["VertexLB"],
            u1=0.0,
            u2=UNSET,
            ur3=UNSET,
        )
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="RigidConsE112",
            region=self.assembly.instances[self.instance].sets["VertexLT"],
            u1=0.0,
            u2=0.0,
            ur3=UNSET,
        )

    def rigid_cons_22(self):
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="RigidConsE221",
            region=self.assembly.instances[self.instance].sets["VertexLB"],
            u1=0.0,
            u2=0.0,
            ur3=UNSET,
        )
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="RigidConsE222",
            region=self.assembly.instances[self.instance].sets["VertexRB"],
            u1=UNSET,
            u2=0.0,
            ur3=UNSET,
        )

    def create_path_load(self, load_magnitude, path_table, time_period):
        """create path loads

        Parameters
        ----------
        load_magnitude : list
            loads magitude, such as [0.02, 0.02,0.02]
        path_table : list
            the exact path
        time_period : float
            simulation time
        """
        # loads should be a list with three elemnts
        E11 = load_magnitude[0]
        E22 = load_magnitude[1]
        E12 = load_magnitude[2]
        num_point = path_table.shape[0]
        # create the path amplitude first
        path_name = ["xx", "yy", "xy"]
        for ii, name in enumerate(path_name):
            path_table_temp = numpy.zeros((num_point, 2))
            path_table_temp[:, 0] = numpy.linspace(
                0, time_period, num_point, endpoint=True
            )
            path_table_temp[:, 1] = path_table[:, ii]
            self._create_random_path(path_table_temp, name)

        # create load for corresonding part
        # create E_11
        self.model.DisplacementBC(
            name="E11",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=E11,
            u2=UNSET,
            ur3=UNSET,
            amplitude=path_name[0],
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
        # create E_22
        self.model.DisplacementBC(
            name="E22",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-T"],
            u1=UNSET,
            u2=E22,
            ur3=UNSET,
            amplitude=path_name[1],
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
        # create E_12
        self.model.DisplacementBC(
            name="E12",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=UNSET,
            u2=E12,
            ur3=UNSET,
            amplitude=path_name[2],
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
        # create E_21
        self.model.DisplacementBC(
            name="E21",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-T"],
            u1=E12,
            u2=UNSET,
            ur3=UNSET,
            amplitude=path_name[2],
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )

    def _create_random_path(self, path_tale, path_name):

        # generate the table
        self.model.TabularAmplitude(
            name=path_name,
            timeSpan=TOTAL,
            smooth=SOLVER_DEFAULT,
            data=(path_tale),
        )
