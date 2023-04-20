# abaqus
# import python abaqus libraries

from abaqus import *
from abaqusConstants import *
from caeModules import *


class CommonProcedure(object):
    def delete_existed_models(self):
        """
        This function is used to delete the existed abaqus models
        :param None
        :return: None
        """
        if "Model-1" in mdb.models.keys():
            del mdb.models["Model-1"]

    def create_new_model(self):
        """
        create a new abaqus model
        :param case_name: a tring of the new model
        :return: None
        """
        Mdb()
        self.model = mdb.Model(
            name=self.model_name, modelType=STANDARD_EXPLICIT
        )

    def delta_for_mesh(self):
        """
        calculate the
        param: none
        return: delta, a float to loose the selection of edges and vertices
        """
        delta = min(min(self.length, self.width) / 1000, self.mesh_size / 10)
        return delta


    def create_sets_for_edges(self):
        """
        an internal function used to create the edges sets of 2D RVE

        param: none
        return: the names of the edges of the RVE

        """
        # information for delta
        delta = self.delta_for_mesh()
        # create sets for edges
        s = self.part.edges
        edgesLEFT = s.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] - self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        self.part.Set(edges=edgesLEFT, name="edgesLEFT")
        edgesRIGHT = s.getByBoundingBox(
            self.center[0] + self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        self.part.Set(edges=edgesRIGHT, name="edgesRIGHT")
        edgesTOP = s.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] + self.width / 2 - delta,
            0,
            self.center[0] + self.width / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        self.part.Set(edges=edgesTOP, name="edgesTOP")
        edgesBOT = s.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] - self.width / 2 + delta,
            0,
        )
        self.part.Set(edges=edgesBOT, name="edgesBOT")
        name_edges = ["edgesLEFT", "edgesRIGHT", "edgesTOP", "edgesBOT"]

        return name_edges

    def create_sets_for_vertices(self):
        """
        this function is used to create sets for vertices
        param: none
        return: names of the vertices
        """
        # information for delta
        delta = self.delta_for_mesh()
        v = self.part.vertices
        vertexLB = v.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] - self.length / 2 + delta,
            self.center[1] - self.width / 2 + delta,
            0,
        )
        self.part.Set(vertices=vertexLB, name="VertexLB")
        vertexRB = v.getByBoundingBox(
            self.center[0] + self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] - self.width / 2 + delta,
            0,
        )
        self.part.Set(vertices=vertexRB, name="VertexRB")
        vertexRT = v.getByBoundingBox(
            self.center[0] + self.length / 2 - delta,
            self.center[1] + self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        self.part.Set(vertices=vertexRT, name="VertexRT")
        vertexLT = v.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] + self.width / 2 - delta,
            0,
            self.center[0] - self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        self.part.Set(vertices=vertexLT, name="VertexLT")
        vertices_name = ["VertexLB", "VertexLT", "VertexRB", "VertexRT"]

        return vertices_name

    def create_assembly(self):

        if "Final_Stuff-1" in self.model.rootAssembly.features.keys():
            self.model.rootAssembly.features.changeKey(
                fromName="Final_Stuff-1", toName=self.instance_name
            )
        else:
            self.assembly = self.model.rootAssembly

        self.assembly.Instance(
            name=self.instance_name, part=self.part, dependent=ON
        )
        self.part = self.model.parts[self.part_name]
