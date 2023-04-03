#                                                                       Modules
# =============================================================================
# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *
## local
from jobs import AbaqusJobs
from loadings import Loading2D
from odbAccess import *
from pbc import PBC2D

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================

# =============================================================================


class RVE2DBase(object):
    def delete_existed_models(self):
        """delete the existed model"""

        if "Model-1" in mdb.models.keys():
            del mdb.models["Model-1"]

    def create_new_model(self):
        """create a new abaqus model"""
        Mdb()
        self.model = mdb.Model(
            name=self.model_name, modelType=STANDARD_EXPLICIT
        )

    def delta_for_mesh(self):
        """calculate the minimum length for mesh"""
        delta = min(min(self.length, self.width) / 1000, self.mesh_size / 10)
        return delta

    def create_sets_for_edges(self):
        """create set for edges"""
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
        """create set for verices"""
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
        name_vertex = ["VertexLB", "VertexLT", "VertexRB", "VertexRT"]

        return name_vertex

    def create_assembly(self):
        """create an assembly for the abaqus model"""
        if "Final_Stuff-1" in self.model.rootAssembly.features.keys():
            self.model.rootAssembly.features.changeKey(
                fromName="Final_Stuff-1", toName="Final_Stuff"
            )
            self.assembly = self.model.rootAssembly
            self.part = self.model.parts["Final_Stuff"]
        else:
            self.assembly = self.model.rootAssembly
            self.assembly.Instance(
                name=self.instance_name, part=self.part, dependent=ON
            )

    def creat_coarse_mesh(self, name_edges):
        """
            mesh the 2D RVE model. In this function, the edges of 2D RVE will be
            divided into certain portions. However, the edges of the fibers
            (inclusions) will be simply given similar size.
            The plane strain element is used.

        Parameters
        ----------
        name_edges : list
            name of edges
        """

        import mesh

        elemType1 = mesh.ElemType(
            elemCode=CPE8R,
            elemLibrary=STANDARD,
            secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED,
            distortionControl=DEFAULT,
        )
        elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
        self.part.setElementType(
            regions=(self.part.faces[:],), elemTypes=(elemType1, elemType2)
        )
        self.part.seedPart(
            size=self.mesh_size, deviationFactor=0.4, minSizeFactor=0.4
        )
        for ii in range(len(name_edges)):
            self.part.seedEdgeBySize(
                edges=self.part.sets[name_edges[ii]].edges,
                size=self.mesh_size,
                deviationFactor=0.4,
                constraint=FIXED,
            )
        self.part.generateMesh()

    def create_2d_pbc(self):
        """
        Assign the periodical boundary condition to the 2D RVE. Since the pbc
        is general for most of the 2D RVEs,this function can be used for most
        of the 2D RVE simulation.

        """

        pbc_generator = PBC2D(
            self.model, self.part, self.assembly, self.length, self.width
        )
        pbc_generator.create_ref_point(self.length, self.width, self.center)
        pbc_generator.vertex_pbc(
            ["NodeLB", "NodeLT", "NodeRB", "NodeRT"],
            ["VertexLB", "VertexLT", "VertexRB", "VertexRT"],
            name_instance=self.instance_name,
        )
        pbc_generator.edges_pbc(name_instance=self.instance_name)

    def create_loads(self, loads):
        """create loads with a constant value

        Parameters
        ----------
        loads : list
            loads 
        """

        loading_conditon = Loading2D(
            self.model, self.assembly, self.instance_name
        )
        loading_conditon.create_loads(loads=loads)

    def create_sequential_jobs(self, subroutine_path):
        """create sequential jobs

        Parameters
        ----------
        platform : str
            simulation platform
        subroutine_path : str
            path of user subroutine, by default ''
        """
        JobGenerator = AbaqusJobs(
            model_name=self.model_name,
            job_name=self.job_name,
            num_cpu=self.num_cpu,
            subtoutine_path=subroutine_path,
        )
        JobGenerator.sequential_job()
