# extend the system path
import sys

# import python libraries
import numpy
# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *
# import packages for abaqus post-processing
from odbAccess import *

try:
    import cPickle as pickle
except ValueError:
    import pickle

from base import RVE2DBase
from geometry import HeterHollowPlate
from loadings import Loading2D
from material import AbaqusMaterialLib
from postprocess import RVEPostProcess2D


class TaylanRVE(RVE2DBase):
    def __init__(
        self,
        sim_info={
            "length": 1.0,
            "width": 1.0,
            "radius": 0.1,
            "location_information": None,
            "youngs_modulus": 200000,
            "poission_ratio": 0.3,
            "yield_table": None,
            "mesh_partition": 100,
            "loads": [0.02, 0.02, 0.02],
            "loads_path": None,
            "time_period": 1.0,
            "job_name": "taylanrve",
            "num_cpu": 1,
            "platform": "ubuntu",
        },
    ):
        """
        Initializaton function of Hollow Plate RVE case

        """

        # define the names of RVE simulation and set them as private variables
        self.model_name = "RVE"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = str(sim_info["job_name"])
        self.platform = sim_info["platform"]
        self.num_cpu = sim_info["num_cpu"]

        # define the import elements of RVE
        self.model = None
        self.sketch = None
        self.part = None
        self.assembly = None
        self.material = None
        #  define the geometry information
        self.length = sim_info["length"]
        self.width = sim_info["width"]
        self.location_information = sim_info["location_information"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )
        self.center = [
            (sim_info["length"]) / 2.0,
            (sim_info["width"]) / 2.0,
        ]
        # loading condition
        self.loads = sim_info["loads"]

        self.youngs_modulus = sim_info["youngs_modulus"]
        self.poission_ratio = sim_info["poission_ratio"]
        self.time_period = sim_info["time_period"]

        # opreation on transfer list to numpy
        num_yield_points = len(sim_info["yield_table"][0])
        self.yield_criterion = numpy.zeros((num_yield_points, 2))
        for ii in range(2):
            self.yield_criterion[:, ii] = sim_info["yield_table"][ii]

        num_path_points = len(sim_info["loads_path"][0])
        self.loads_path = numpy.zeros((num_path_points, 3))
        for ii in range(3):
            self.loads_path[:, ii] = sim_info["loads_path"][ii]

        # run the simulation
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
        self._create_step(timePeriod=self.time_period)
        self._create_path_load()
        self.create_sequential_jobs(subroutine_path="")
        if self.platform == "cluster":
            RVEPostProcess2D(self.job_name)

    def _create_part(self):
        PartGenerator = HeterHollowPlate(
            length=self.length,
            width=self.width,
            location_info=self.location_information,
            name_part=self.part_name,
            model=self.model,
        )
        self.part = PartGenerator.create_part()

    def _create_material(self, name_faces):
        # use plasticity material
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Matrix",
            model=self.model,
            part=self.part,
            name_set=name_faces[0],
        )
        MaterialGenerator.CreateVonMisesPlasticMaterial(
            E=self.youngs_modulus,
            v=self.poission_ratio,
            yield_criterion=self.yield_criterion,
        )

    def _create_path_load(self):
        load_creater = Loading2D(self.model, self.assembly, self.instance_name)
        load_creater.create_path_load(
            self.loads, self.loads_path, self.time_period
        )

    def _create_geometry_set(self):

        #  create part for faces
        p = self.model.parts[self.part_name]
        faces = p.faces[:]
        p.Set(faces=faces, name="all_faces")
        name_faces = ["all_faces"]
        # create sets for edges:
        name_edges = self.create_sets_for_edges()
        name_vertex = self.create_sets_for_vertices()
        return name_faces, name_edges, name_vertex

    def _create_step(
        self,
        initialInc=0.01,
        maxInc=1,
        maxNumInc=1000000,
        minInc=1e-20,
        timePeriod=10,
    ):

        self.model.StaticStep(name="Step-1", previous="Initial", nlgeom=ON)
        self.model.StaticStep(
            initialInc=initialInc,
            maxInc=maxInc,
            maxNumInc=maxNumInc,
            minInc=minInc,
            name="Step-1",
            previous="Initial",
            timePeriod=timePeriod,
        )

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
            timeInterval=0.01,
        )
        self.model.FieldOutputRequest(
            name="F-Output-2",
            createStepName="Step-1",
            variables=("U", "RF"),
            timeInterval=0.01,
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
            ),
            timeInterval=0.01,
        )
