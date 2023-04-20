# abaqus
import json

import numpy
from abaqus import *
from abaqusConstants import *
from caeModules import *


class Loading2D:
    pass


class LoadingRigid2D:
    def rigid_constraint_for_complex_loading(self):
        """displacement rigid constraints for complex loading where
        at least shear and tension/compression loadings are applied
        """
        # select two nodes to restrict the rigid movement
        allnodes = self.part.nodes
        xx_min = self.center[0] - self.length / 2
        xx_max = self.center[0] + self.length / 2
        yy_min = self.center[1] - self.width / 2
        yy_max = self.center[1] + self.width / 2
        ysupportnodes = allnodes.getByBoundingBox(
            (0.1 * xx_max + 1.9 * xx_min) / 2,
            (0.1 * yy_max + 1.9 * yy_min) / 2,
            0,
            (0.3 * xx_max + 1.7 * xx_min) / 2,
            (0.3 * yy_max + 1.7 * yy_min) / 2,
            0,
        )
        # assembly set for y support node
        self.assembly.SetFromNodeLabels(
            name="ysupportnode",
            nodeLabels=((self.instance_name, (ysupportnodes[0].label,)),),
            unsorted=True,
        )
        xsupportnodes = allnodes.getByBoundingBox(
            (1.7 * xx_max + 0.3 * xx_min) / 2,
            (1.7 * yy_max + 0.3 * yy_min) / 2,
            0,
            (1.9 * xx_max + 0.1 * xx_min) / 2,
            (1.9 * yy_max + 0.1 * yy_min) / 2,
            0,
        )
        # assembly set for x support node
        self.assembly.SetFromNodeLabels(
            name="xsupportnode",
            nodeLabels=((self.instance_name, (xsupportnodes[0].label,)),),
            unsorted=True,
        )
        # rigid movement restriction
        # x direction
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="rigid_x",
            region=self.assembly.sets["xsupportnode"],
            u1=0.0,
            u2=UNSET,
            ur3=UNSET,
        )
        # x direction
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="rigid_y",
            region=self.assembly.sets["ysupportnode"],
            u1=UNSET,
            u2=0.0,
            ur3=UNSET,
        )

    def rigid_constraints_for_x_direction(self):
        """displacement rigid constraints for x direction"""
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="rigid_x_1",
            region=self.assembly.instances[self.instance_name].sets["VertexLB"],
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
            name="rigid_x_2",
            region=self.assembly.instances[self.instance_name].sets["VertexLT"],
            u1=0.0,
            u2=0.0,
            ur3=UNSET,
        )

    def rigid_constraints_for_y_direction(self):
        self.model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="rigid_y_1",
            region=self.assembly.instances[self.instance_name].sets["VertexLB"],
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
            name="rigid_y_2",
            region=self.assembly.instances[self.instance_name].sets["VertexRB"],
            u1=UNSET,
            u2=0.0,
            ur3=UNSET,
        )


class NormalDisplacementLoading(LoadingRigid2D):

    def create_load(self):
        E11 = self.strain[0]
        E22 = self.strain[1]
        E12 = self.strain[2]

        if E22 == E12 == 0.0:
            # x direction tension/compression
            self.E11(E11=E11)
            # rigid constrains
            self.rigid_constraints_for_x_direction()
        elif E11 == E12 == 0.0:
            # y direction tension/compression
            self.E22(E22=E22)
            # rigid movement
            self.rigid_constraints_for_y_direction()
        elif E11 == E22 == 0.0:
            # shear loading
            self.E12(E12=E12)
            # rigid movement
            self.rigid_constraints_for_x_direction()
        else:
            # complex loading
            self.E11(E11=E11)
            self.E12(E12=E12)
            self.E22(E22=E22)
            # rigid movement
            self.rigid_constraint_for_complex_loading()

    def E11(self, E11):
        "displacement for x direction"
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

    def E22(self, E22):
        "displacement for y direction"
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

    def E12(self, E12):
        "displacement for xy direction"
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


class HistoryDependentDisplacement2D(LoadingRigid2D):
 

    def create_load(self):

        # loads should be a list with three elements
        E11 = self.strain[0]
        E22 = self.strain[1]
        E12 = self.strain[2]
        num_point = self.strain_amplitude.shape[0]
        # create the path amplitude first
        path_name = ["xx", "yy", "xy"]
        for ii, name in enumerate(path_name):
            path_table_temp = numpy.zeros((num_point, 2))
            path_table_temp[:, 0] = numpy.linspace(
                0, self.time_period, num_point, endpoint=True
            )
            path_table_temp[:, 1] = self.strain_amplitude[:, ii]
            self._create_random_path(path_table_temp, name)

        if E22 == E12 == 0.0:
            # x direction tension/compression
            self.E11(E11=E11, path_name=path_name[0])
            # rigid constrains
            self.rigid_constraints_for_x_direction()
        elif E11 == E12 == 0.0:
            # y direction tension/compression
            self.E22(E22=E22, path_name=path_name[1])
            # rigid movement
            self.rigid_constraints_for_y_direction()
        elif E11 == E22 == 0.0:
            # shear loading
            self.E12(E12=E12, path_name=path_name[2])
            # rigid movement
            self.rigid_constraints_for_x_direction()
        else:
            # complex loading
            self.E11(E11=E11, path_name=path_name[0])
            self.E22(E22=E22, path_name=path_name[1])
            self.E12(E12=E12, path_name=path_name[2])
            # rigid movement
            self.rigid_constraint_for_complex_loading()

    def E11(self,E11, path_name):
        # create load for corresponding part
        # create E_11
        self.model.DisplacementBC(
            name="E11",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=E11,
            u2=UNSET,
            ur3=UNSET,
            amplitude=path_name,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
    def E22(self, E22, path_name):
        # create E_22
        self.model.DisplacementBC(
            name="E22",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-T"],
            u1=UNSET,
            u2=E22,
            ur3=UNSET,
            amplitude=path_name,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
    def E12(self,E12, path_name): 
        # create E_12
        self.model.DisplacementBC(
            name="E12",
            createStepName="Step-1",
            region=self.assembly.sets["Ref-R"],
            u1=UNSET,
            u2=E12,
            ur3=UNSET,
            amplitude=path_name,
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
            amplitude=path_name,
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