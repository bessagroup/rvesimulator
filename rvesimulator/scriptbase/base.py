# abaqus
# import python libraries
import numpy
from abaqus import *
from abaqusConstants import *
from caeModules import *
## import packages for abaqus post-processing
from odbAccess import *

try:
    import cPickle as pickle  # Improve speed
except ValueError:
    import pickle

from jobs import AbaqusJobs
from loadings import Loading2D
from pbc import PBC2D


class RVE2DBase(object):

    def delete_existed_models(self):
        """
        This function is used to delete the existed abaqus models 
        :param None
        :return: None
        """
        if 'Model-1' in mdb.models.keys():
            del mdb.models['Model-1']

    def create_new_model(self):
        """
        create a new abaqus model 
        :param case_name: a tring of the new model
        :return: None 
        """
        Mdb()
        self.model = mdb.Model(name=self.model_name, modelType=STANDARD_EXPLICIT)

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
        edgesLEFT = s.getByBoundingBox(self.center[0] - self.length / 2 - delta,
                                       self.center[1] - self.width / 2 - delta, 0,
                                       self.center[0] - self.length / 2 + delta,
                                       self.center[1] + self.width / 2 + delta, 0)
        self.part.Set(edges=edgesLEFT, name='edgesLEFT')
        edgesRIGHT = s.getByBoundingBox(self.center[0] + self.length / 2 - delta,
                                        self.center[1] - self.width / 2 - delta, 0,
                                        self.center[0] + self.length / 2 + delta,
                                        self.center[1] + self.width / 2 + delta, 0)
        self.part.Set(edges=edgesRIGHT, name='edgesRIGHT')
        edgesTOP = s.getByBoundingBox(self.center[0] - self.length / 2 - delta, self.center[1] + self.width / 2 - delta,
                                      0,
                                      self.center[0] + self.width / 2 + delta, self.center[1] + self.width / 2 + delta,
                                      0)
        self.part.Set(edges=edgesTOP, name='edgesTOP')
        edgesBOT = s.getByBoundingBox(self.center[0] - self.length / 2 - delta, self.center[1] - self.width / 2 - delta,
                                      0,
                                      self.center[0] + self.length / 2 + delta, self.center[1] - self.width / 2 + delta,
                                      0)
        self.part.Set(edges=edgesBOT, name='edgesBOT')
        name_edges = ['edgesLEFT', 'edgesRIGHT', 'edgesTOP', 'edgesBOT']

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
        vertexLB = v.getByBoundingBox(self.center[0] - self.length / 2 - delta, self.center[1] - self.width / 2 - delta,
                                      0,
                                      self.center[0] - self.length / 2 + delta, self.center[1] - self.width / 2 + delta,
                                      0)
        self.part.Set(vertices=vertexLB, name='VertexLB')
        vertexRB = v.getByBoundingBox(self.center[0] + self.length / 2 - delta, self.center[1] - self.width / 2 - delta,
                                      0,
                                      self.center[0] + self.length / 2 + delta, self.center[1] - self.width / 2 + delta,
                                      0)
        self.part.Set(vertices=vertexRB, name='VertexRB')
        vertexRT = v.getByBoundingBox(self.center[0] + self.length / 2 - delta, self.center[1] + self.width / 2 - delta,
                                      0,
                                      self.center[0] + self.length / 2 + delta, self.center[1] + self.width / 2 + delta,
                                      0)
        self.part.Set(vertices=vertexRT, name='VertexRT')
        vertexLT = v.getByBoundingBox(self.center[0] - self.length / 2 - delta, self.center[1] + self.width / 2 - delta,
                                      0,
                                      self.center[0] - self.length / 2 + delta, self.center[1] + self.width / 2 + delta,
                                      0)
        self.part.Set(vertices=vertexLT, name='VertexLT')
        name_vertex = ['VertexLB', 'VertexLT', 'VertexRB', 'VertexRT']

    def create_assembly(self):
        """
        create an assembly for the abaqus model 
        param: none
        return: none
        """
        if 'Final_Stuff-1' in self.model.rootAssembly.features.keys():
            self.model.rootAssembly.features.changeKey(fromName='Final_Stuff-1', toName='Final_Stuff')
            self.assembly = self.model.rootAssembly
            self.part = self.model.parts['Final_Stuff']
        else:
            self.assembly = self.model.rootAssembly
            self.assembly.Instance(name=self.instance_name, part=self.part, dependent=ON)

    def creat_coarse_mesh(self, name_edges):
        """
        mesh the 2D RVE model. In this function, the edges of 2D RVE will be divided into 
        certain portions. However, the edges of the fibers(inclusions) will be simply given similar size. 
        The plane strain element is used. 

        param: name_edges, a list that contains the names of the edges.  

        return: none 
        """

        import mesh
        elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD,
                                  secondOrderAccuracy=OFF, hourglassControl=ENHANCED,
                                  distortionControl=DEFAULT)
        elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
        self.part.setElementType(regions=(self.part.faces[:],), elemTypes=(elemType1, elemType2))
        self.part.seedPart(size=self.mesh_size, deviationFactor=0.4, minSizeFactor=0.4)
        for ii in range(len(name_edges)):
            self.part.seedEdgeBySize(edges=self.part.sets[name_edges[ii]].edges, size=self.mesh_size,
                                     deviationFactor=0.4, constraint=FIXED)
        self.part.generateMesh()

    def create_2d_pbc(self):

        """
        Assign the periodical boundary condition to the 2D RVE. Since the pbc is general for most 
        of the 2D RVEs,this function can be used for most of the 2D RVE simulation.  
        Note: 

        param: none
        return: none 

        """

        pbc_generator = PBC2D(self.model, self.part, self.assembly, self.length, self.width)
        pbc_generator.create_ref_point(self.length, self.width, self.center)
        pbc_generator.vertex_pbc(['NodeLB', 'NodeLT', 'NodeRB', 'NodeRT'],
                                 ['VertexLB', 'VertexLT', 'VertexRB', 'VertexRT'], name_instance=self.instance_name)
        pbc_generator.edges_pbc(name_instance=self.instance_name)

    def create_loads(self, loads):
        """This functions is used to create the loads for 2D RVE 

        Args:
            loads (list): a list contains the loads, defualt: loads=[0.5, 0.0, 0.0]
        Return: 
            None
        """
        loading_conditon = Loading2D(self.model, self.assembly, self.instance_name)
        loading_conditon.create_loads(loads=loads)

    def create_sequential_jobs(self, subroutine_path):
        JobGenerator = AbaqusJobs(model_name=self.model_name,
                                  job_name=self.job_name,
                                  subtoutine_path=subroutine_path)
        JobGenerator.sequential_job()
