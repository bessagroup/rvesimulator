# abaqus
import json

import numpy
from abaqus import *
from abaqusConstants import *
from caeModules import *


class Mesh2D:
    def create_mesh(self, edges_name, element_type):
        import mesh

        # get the element type
        if element_type == "linear":
            elemType1, elemType2 = self._define_linear_element()
        elif element_type == "quadratic":
            elemType1, elemType2 = self._define_quadratic_element()

        # set all part with defined element type
        self.part.setElementType(
            regions=(self.part.faces[:],), elemTypes=(elemType1, elemType2)
        )
        # overall seed (for the circles)
        self.part.seedPart(
            size=self.mesh_size, deviationFactor=0.4, minSizeFactor=0.4
        )

        # mesh every edge (for all edges)
        for ii in range(len(edges_name)):
            self.part.seedEdgeBySize(
                edges=self.part.sets[edges_name[ii]].edges,
                size=self.mesh_size,
                deviationFactor=0.4,
                constraint=FIXED,
            )
        # generate mesh
        self.part.generateMesh()

    def _define_linear_element(self):

        elemType1 = mesh.ElemType(
            elemCode=CPE4R,
            elemLibrary=STANDARD,
            secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED,
            distortionControl=DEFAULT,
        )
        elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)

        return elemType1, elemType2

    def _define_quadratic_element(self):

        elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)

        return elemType1, elemType2
