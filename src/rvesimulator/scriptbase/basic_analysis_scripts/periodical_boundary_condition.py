# abaqus
# third-party
import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *


class PeriodicalBoundaryCondition2D:
    def create_reference_points(self):
        # right reference point
        right_reference_point_id = self.assembly.ReferencePoint(
            point=(self.length * 1.5, self.center[1], 0.0)
        ).id
        # top reference point
        top_reference_point_id = self.assembly.ReferencePoint(
            point=(self.center[0], self.width * 1.5, 0.0)
        ).id
        # reference points
        reference_points = self.assembly.referencePoints
        # create sets for reference points
        self.assembly.Set(
            name="Ref-R",
            referencePoints=((reference_points[right_reference_point_id],)),
        )
        self.assembly.Set(
            name="Ref-T",
            referencePoints=((reference_points[top_reference_point_id],)),
        )

    def pbc_for_vertices(self, set_name, geometry_name):

        # import session
        import assembly

        # regenerate assembly
        session.viewports["Viewport: 1"].setValues(
            displayedObject=self.assembly
        )
        self.assembly.regenerate()

        for ii in range(len(geometry_name)):
            vertex_node = self.part.sets[geometry_name[ii]].nodes
            # create the assembly node set for vertex
            self.assembly.SetFromNodeLabels(
                name=set_name[ii],
                nodeLabels=((self.instance_name, (vertex_node[0].label,)),),
                unsorted=True,
            )

        # apply the pbc to vertexes
        self.model.Equation(
            name="LB_RT_1",
            terms=(
                (1, "NodeRT", 1),
                (-1, "NodeLB", 1),
                (-1 * self.length, "Ref-R", 1),
                (-1 * self.width, "Ref-T", 1),
            ),
        )
        self.model.Equation(
            name="LB_RT_2",
            terms=(
                (1, "NodeRT", 2),
                (-1, "NodeLB", 2),
                (-1 * self.length, "Ref-R", 2),
                (-1 * self.width, "Ref-T", 2),
            ),
        )
        self.model.Equation(
            name="LT_RB_1",
            terms=(
                (1, "NodeRB", 1),
                (-1, "NodeLT", 1),
                (-1 * self.length, "Ref-R", 1),
                (1 * self.width, "Ref-T", 1),
            ),
        )
        self.model.Equation(
            name="LT_RB_2",
            terms=(
                (1, "NodeRB", 2),
                (-1, "NodeLT", 2),
                (-1 * self.length, "Ref-R", 2),
                (1 * self.width, "Ref-T", 2),
            ),
        )

    def pbc_for_edges(self):
        import assembly

        session.viewports["Viewport: 1"].setValues(
            displayedObject=self.assembly
        )
        self.assembly.regenerate()
        # part 1: equations for edges 2 (edgesFRONT_RIGHT) and 4 (edgesBACK_LEFT)
        edgesRIGHT_nodes = self.part.sets["edgesRIGHT"].nodes
        edgesRIGHT_nodes_sorted = sorted(edgesRIGHT_nodes, key=get_node_y)
        edgesLEFT_nodes = self.part.sets["edgesLEFT"].nodes
        edgesLEFT_nodes_sorted = sorted(edgesLEFT_nodes, key=get_node_y)
        if len(edgesRIGHT_nodes_sorted) == len(edgesLEFT_nodes_sorted):
            for ii in range(1, len(edgesRIGHT_nodes_sorted) - 1):
                self.assembly.SetFromNodeLabels(
                    name="LEFT_" + str(ii),
                    nodeLabels=(
                        (
                            self.instance_name,
                            tuple([edgesLEFT_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                self.assembly.SetFromNodeLabels(
                    name="RIGHT_" + str(ii),
                    nodeLabels=(
                        (
                            self.instance_name,
                            tuple([edgesRIGHT_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                for jj in range(1, 3):
                    self.model.Equation(
                        name="LEFT_RIGHT_" + str(ii) + "_" + str(jj),
                        terms=(
                            (1, "RIGHT_" + str(ii), jj),
                            (-1, "LEFT_" + str(ii), jj),
                            (-1 * self.length, "Ref-R", jj),
                        ),
                    )
        else:
            print
            "the number of nodes between the two sides are not the same"
        # part II:
        edgesTOP_nodes = self.part.sets["edgesTOP"].nodes
        edgesTOP_nodes_sorted = sorted(edgesTOP_nodes, key=get_node_x)
        edgesBOT_nodes = self.part.sets["edgesBOT"].nodes
        edgesBOT_nodes_sorted = sorted(edgesBOT_nodes, key=get_node_x)
        if len(edgesTOP_nodes_sorted) == len(edgesBOT_nodes_sorted):
            for ii in range(1, len(edgesTOP_nodes_sorted) - 1):
                self.assembly.SetFromNodeLabels(
                    name="TOP_" + str(ii),
                    nodeLabels=(
                        (
                            self.instance_name,
                            tuple([edgesTOP_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                self.assembly.SetFromNodeLabels(
                    name="BOT_" + str(ii),
                    nodeLabels=(
                        (
                            self.instance_name,
                            tuple([edgesBOT_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                for jj in range(1, 3):
                    self.model.Equation(
                        name="TOP_BOT_" + str(ii) + "_" + str(jj),
                        terms=(
                            (1, "TOP_" + str(ii), jj),
                            (-1, "BOT_" + str(ii), jj),
                            (-1 * self.width, "Ref-T", jj),
                        ),
                    )
        else:
            print
            "the number of nodes between the two sides are not the same"


def get_node_y(node):
    return node.coordinates[1]


def get_node_x(node):
    return node.coordinates[0]
