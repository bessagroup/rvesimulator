# extend the system path
import sys

# import python libraries
import numpy

# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *

## import packages for abaqus post-processing
from odbAccess import *

try:
    import cPickle as pickle  # Improve speed
except ValueError:
    import pickle

from base import RVE2DBase
from geometry import HollowPlate
from material import AbaqusMaterialLib


class HollowPlateRVE(RVE2DBase):
    def __init__(self, sim_info=None):
        """
        Initializaton function of Hollow Plate RVE case

        """

        # define the names of RVE simulation and set them as private variables
        self.model_name = "RVE"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = str(sim_info["job_name"])
        # define the import elements of RVE
        self.model = None
        self.sketch = None
        self.part = None
        self.assembly = None
        self.material = None
        #  define the geometry information
        self.length = sim_info["size"]
        self.width = sim_info["size"]
        self.center = [0.0, 0.0]
        self.radius = sim_info["radius"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_portion"]
        )

        # loading condition
        self.loads = sim_info["loads"]
        self.youngs_modulus = sim_info["youngs_modulus"]
        self.poission_ratio = sim_info["poission_ratio"]

        # excute the simulation process
        self.ScriptGenerator()

    def ScriptGenerator(self):

        self.create_new_model()
        self.delete_existed_models()
        self._create_part()
        self.create_assembly()
        name_faces, name_edges, name_vertex = self._create_geometry_set()
        self._create_material(name_faces)
        self.creat_coarse_mesh(name_edges)
        self.create_2d_pbc()
        self._create_step()
        self.create_loads(loads=self.loads)
        self.create_sequential_jobs(subroutine_path="")

    def _create_part(self):
        PartGenerator = HollowPlate(
            point1=(-self.length / 2, -self.width / 2),
            point2=(self.length / 2, self.width / 2),
            center=self.center,
            radius=self.radius,
            model=self.model,
            name_part=self.part_name,
        )
        self.part = PartGenerator.create_part()

    def _create_material(self, name_faces):
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Matrix",
            model=self.model,
            part=self.part,
            name_set=name_faces[0],
        )
        self.material = MaterialGenerator.CreateElasticMaterial(
            E=self.youngs_modulus, v=self.poission_ratio
        )

    def _create_geometry_set(self):

        # create part for faces
        p = self.model.parts[self.part_name]
        faces = p.faces[:]
        p.Set(faces=faces, name="all_faces")
        name_faces = ["all_faces"]
        # create sets for edges
        name_edges = self.create_sets_for_edges()
        name_vertex = self.create_sets_for_vertices()

        return name_faces, name_edges, name_vertex

    def _create_step(
        self,
        initialInc=0.1,
        maxInc=1,
        maxNumInc=1000000,
        minInc=1e-20,
        timePeriod=1,
    ):

        self.model.StaticStep(name="Step-1", previous="Initial")
        self.model.StaticStep(
            initialInc=initialInc,
            maxInc=maxInc,
            maxNumInc=maxNumInc,
            minInc=minInc,
            name="Step-1",
            previous="Initial",
            timePeriod=timePeriod,
        )

        ## create Final-outputs
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
            timeInterval=0.1,
        )
        self.model.FieldOutputRequest(
            name="F-Output-2",
            createStepName="Step-1",
            variables=("U", "RF"),
            timeInterval=0.1,
        )
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
            timeInterval=0.1,
        )
