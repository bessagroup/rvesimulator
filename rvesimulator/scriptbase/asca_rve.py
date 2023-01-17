#                                                                       Modules
# =============================================================================
# standard
import sys

# thrid party
import numpy

try:
    import cPickle as pickle  # Improve speed
except ValueError:
    import pickle
# abaqus
from abaqus import *
from abaqusConstants import *

# local
from base import RVE2DBase
from caeModules import *
from geometry import CircleInclusion
from material import AbaqusMaterialLib
from odbAccess import *

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
class ASCARVE(RVE2DBase):
    def __init__(self, sim_info=None):
        # names of model, part, instance
        self.model_name = "RVE"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = str(sim_info["job_name"])
        self.num_cpu = sim_info["num_cpu"]
        self.platform = sim_info["platform"]
        # define the import elements of RVE
        self.model = None
        self.sketch = None
        self.part = None
        self.assembly = None
        self.material = None
        # information of geometry of RVE
        self.loc_info = sim_info
        # mech and sets information
        self.loc = sim_info["location_information"]
        self.length = (
            sim_info["len_end"] - sim_info["len_start"]
        ) - 2 * sim_info["radius_mu"]
        self.width = (
            sim_info["wid_end"] - sim_info["wid_start"]
        ) - 2 * sim_info["radius_mu"]
        self.center = [
            (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
            (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
        ]
        self.radius_mu = sim_info["radius_mu"]
        self.radius_std = sim_info["radius_std"]
        # information of RVE modeling
        self.loads = sim_info["loads"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )
        self.simulation_time = sim_info["simulation_time"]
        # figure out the loading condition

        self.ScriptGenerator()

    def ScriptGenerator(self):
        self.create_new_model()
        self.delete_existed_models()
        self._create_part()
        self.create_assembly()
        name_faces, name_edges, name_vertex = self._create_geometry_set()
        self.creat_coarse_mesh(name_edges)
        self._create_material(name_faces)
        self.create_2d_pbc()
        self._create_step(timePeriod=self.simulation_time)
        self.create_loads(loads=self.loads)
        self.create_sequential_jobs(subroutine_path="")

    def _create_part(self):
        PartGenerator = CircleInclusion(
            self.model,
            self.loc_info,
            self.model_name,
            self.part_name,
            self.instance_name,
        )
        PartGenerator.create_part()

    def _create_geometry_set(self):

        # create part set for faces
        p = self.model.parts[self.part_name]
        faces = p.faces[:]
        p.Set(faces=faces, name="all_faces")
        fiberface = p.faces.getByBoundingCylinder(
            (self.loc[0][0], self.loc[0][1], 0.0),
            (self.loc[0][0], self.loc[0][1], 1.0),
            self.loc[0][2] + 0.001 * self.loc[0][2],
        )
        p.Set(faces=fiberface, name="fiberface")
        for ii in range(1, len(self.loc)):
            fiberface_1 = p.faces.getByBoundingCylinder(
                (self.loc[ii][0], self.loc[ii][1], 0.0),
                (self.loc[ii][0], self.loc[ii][1], 1.0),
                self.loc[ii][2] + 0.001 * self.loc[ii][2],
            )
            p.Set(faces=fiberface_1, name="fiberface_1")
            p.SetByBoolean(
                name="fiberface",
                sets=(p.sets["fiberface_1"], p.sets["fiberface"]),
                operation=UNION,
            )
        # delete fiberface_1
        del self.part.sets["fiberface_1"]
        p.SetByBoolean(
            name="matrixface",
            sets=(p.sets["all_faces"], p.sets["fiberface"]),
            operation=DIFFERENCE,
        )
        name_faces = ["matrixface", "fiberface"]

        # create sets for edges:
        name_edges = self.create_sets_for_edges()
        name_vertex = self.create_sets_for_vertices()

        return name_faces, name_edges, name_vertex

    def _create_material(self, name_faces):
        # print name_faces
        # create material for the matrix
        # create the yield criterion for the matrix material
        yield_criterion = numpy.zeros((101, 2))
        yield_criterion[:, 1] = numpy.linspace(0, 1, 101)
        yield_criterion[:, 0] = 0.5 + 0.2 * (yield_criterion[:, 1]) ** 0.4
        yield_criterion[-1, 1] = 10.0
        yield_criterion[-1, 0] = 0.5 + 0.2 * (yield_criterion[-1, 1]) ** 0.4
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Matrix",
            model=self.model,
            part=self.part,
            name_set=name_faces[0],
        )
        # MaterialGenerator.CreateElasticMaterial(E=100.0, v=0.1)
        MaterialGenerator.CreateVonMisesPlasticMaterial(
            E=100, v=0.3, yield_criterion=yield_criterion
        )
        # create the yield criterion for the matrix material
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Fiber",
            model=self.model,
            part=self.part,
            name_set=name_faces[1],
        )
        MaterialGenerator.CreateElasticMaterial(E=1.0, v=0.19)

    def _create_step(
        self,
        initialInc=0.1,
        maxInc=1,
        maxNumInc=1000000,
        minInc=1e-20,
        timePeriod=10,
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
            nlgeom=ON,
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
                "ETOTAL",
            ),
            timeInterval=0.01,
        )
