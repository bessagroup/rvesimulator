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
from geometry import CircleInclusion, FullCircleInclusion, HollowPlate
from loadings import Loading2D
from material import AbaqusMaterialLib
from postprocess import RVEPostProcess2D


class PnasHollowPlate(RVE2DBase):
    def __init__(
        self,
        sim_info={
            "length": 1.0,
            "width": 1.0,
            "radius": 0.2,
            "youngs_modulus": 100.0,
            "poission_ratio": 0.3,
            "yield_table": None,
            "mesh_partition": 30,
            "strain": [0.1, 0.0, 0.0],
            "strain_path": None,
            "time_period": 1,
            "platform": "ubuntu",
            "num_cpu": 1,
        },
    ):
        """
        Initializaton function of Hollow Plate RVE case

        """

        # define the names of RVE simulation and set them as private variables
        self.model_name = "RVE"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = "pnas_hollow_plate"
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
        self.center = [0.0, 0.0]
        self.radius = sim_info["radius"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )

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
        initialInc=0.1,
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


class PnasCompositeRVE(RVE2DBase):
    def __init__(
        self,
        sim_info={
            "location_information": None,
            "len_start": None,
            "len_end": None,
            "wid_start": None,
            "wid_end": None,
            "radius_mu": None,
            "radius_std": None,
            "job_name": "pnas_composite",
            "loads": [0.02, 0.02, 0.02],
            "loads_path": None,
            "E_matrix": None,
            "Pr_matrix": None,
            "yield_table_matrix": None,
            "E_fiber": None,
            "Pr_fiber": None,
            "mesh_partition": None,
            "time_period": 1.0,
            "platform": "ubuntu",
            "num_cpu": 1,
        },
    ):
        # names of model, part, instance
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
        self.time_period = sim_info["time_period"]

        # material properties
        self.E_matrix = sim_info["E_matrix"]
        self.Pr_matrix = sim_info["Pr_matrix"]
        self.E_fiber = sim_info["E_fiber"]
        self.Pr_fiber = sim_info["Pr_fiber"]

        # opreation on transfer list to numpy
        num_yield_points = len(sim_info["yield_table_matrix"][0])
        self.yield_crtierion_matrix = numpy.zeros((num_yield_points, 2))
        for ii in range(2):
            self.yield_crtierion_matrix[:, ii] = sim_info[
                "yield_table_matrix"
            ][ii]

        num_path_points = len(sim_info["loads_path"][0])
        self.loads_path = numpy.zeros((num_path_points, 3))
        for ii in range(3):
            self.loads_path[:, ii] = sim_info["loads_path"][ii]

        self.ScriptGenerator()

    def ScriptGenerator(self):
        self.create_new_model()
        self.delete_existed_models()
        self._create_part()
        self._create_assembly()
        name_faces, name_edges, name_vertex = self._create_geometry_set()
        self._creat_mesh(name_edges)
        self._create_material(name_faces)
        self.create_2d_pbc()
        self._create_step(timePeriod=self.time_period)
        self._create_path_load()
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

    def _create_assembly(self):
        self.model.rootAssembly.features.changeKey(
            fromName="Final_Stuff-1", toName="Final_Stuff"
        )
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Final_Stuff"]

    def _create_geometry_set(self):

        # create part set for faces
        p = self.model.parts[self.part_name]
        faces = p.faces[:]
        p.Set(faces=faces, name="all_faces")
        fiberface = p.faces.getByBoundingCylinder(
            (self.loc[0][0], self.loc[0][1], 0.0),
            (self.loc[0][0], self.loc[0][1], 1.0),
            self.radius_mu + 0.001 * self.radius_mu,
        )
        p.Set(faces=fiberface, name="fiberface")
        for ii in range(1, len(self.loc)):
            fiberface_1 = p.faces.getByBoundingCylinder(
                (self.loc[ii][0], self.loc[ii][1], 0.0),
                (self.loc[ii][0], self.loc[ii][1], 1.0),
                self.radius_mu + 0.001 * self.radius_mu,
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

    def _creat_mesh(self, name_edges):

        import mesh

        elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(
            elemCode=CPE6M,
            elemLibrary=STANDARD,
            secondOrderAccuracy=ON,
            distortionControl=DEFAULT,
        )
        self.part.setElementType(
            regions=(self.part.faces[:],), elemTypes=(elemType1, elemType2)
        )
        self.part.seedPart(
            size=self.mesh_size, deviationFactor=0.2, minSizeFactor=0.9
        )
        for ii in range(len(name_edges)):
            self.part.seedEdgeBySize(
                edges=self.part.sets[name_edges[ii]].edges,
                size=self.mesh_size,
                deviationFactor=0.9,
                constraint=FIXED,
            )
        self.part.generateMesh()

    def _create_material(self, name_faces):
        # assign material property to matrix material
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Matrix",
            model=self.model,
            part=self.part,
            name_set=name_faces[0],
        )
        MaterialGenerator.CreateVonMisesPlasticMaterial(
            E=self.E_matrix,
            v=self.Pr_matrix,
            yield_criterion=self.yield_crtierion_matrix,
        )
        # assign material property to fiber material
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Fiber",
            model=self.model,
            part=self.part,
            name_set=name_faces[1],
        )
        MaterialGenerator.CreateElasticMaterial(
            E=self.E_fiber, v=self.Pr_fiber
        )

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
            nlgeom=ON,
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
                "ETOTAL",
            ),
            timeInterval=0.01,
        )

    def _create_path_load(self):
        load_creater = Loading2D(self.model, self.assembly, self.instance_name)
        load_creater.create_path_load(
            self.loads, self.loads_path, self.time_period
        )


class StatisticRepresentVolume(RVE2DBase):
    def __init__(
        self,
        sim_info={
            "location_information": None,
            "len_start": None,
            "len_end": None,
            "wid_start": None,
            "wid_end": None,
            "radius": None,
            "job_name": "pnas_sve",
            "loads": [0.1, 0.0, 0.0],
            "loads_path": None,
            "E_matrix": None,
            "Pr_matrix": None,
            "yield_table_matrix": None,
            "E_fiber": None,
            "Pr_fiber": None,
            "mesh_partition": 50,
            "time_period": 1.0,
            "platform": "ubuntu",
            "num_cpu": 1,
        },
    ):

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
        ) - 2 * sim_info["radius"]
        self.width = (
            sim_info["wid_end"] - sim_info["wid_start"]
        ) - 2 * sim_info["radius"]
        self.center = [
            (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
            (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
        ]
        self.radius = sim_info["radius"]

        # information of RVE modeling
        self.loads = sim_info["loads"]
        # self.loads_path = sim_info["loads_path"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )
        self.time_period = sim_info["time_period"]

        # material properties
        self.E_matrix = sim_info["E_matrix"]
        self.Pr_matrix = sim_info["Pr_matrix"]
        # self.yield_crtierion_matrix = sim_info["yield_table_matrix"]
        self.E_fiber = sim_info["E_fiber"]
        self.Pr_fiber = sim_info["Pr_fiber"]

        # opreation on transfer list to numpy
        num_yield_points = len(sim_info["yield_table_matrix"][0])
        self.yield_crtierion_matrix = numpy.zeros((num_yield_points, 2))
        for ii in range(2):
            self.yield_crtierion_matrix[:, ii] = sim_info[
                "yield_table_matrix"
            ][ii]

        num_path_points = len(sim_info["loads_path"][0])
        self.loads_path = numpy.zeros((num_path_points, 3))
        for ii in range(3):
            self.loads_path[:, ii] = sim_info["loads_path"][ii]

        self.ScriptGenerator()

    def ScriptGenerator(self):
        self.create_new_model()
        self.delete_existed_models()
        self._create_part()
        self._create_assembly()
        name_faces, name_edges, name_vertex = self._create_geometry_set()
        self._creat_mesh(name_edges)
        self._create_material(name_faces)
        self.create_2d_pbc()
        self._create_step(timePeriod=self.time_period)
        self._create_path_load()
        self.create_sequential_jobs(subroutine_path="")
        if self.platform == "cluster":
            RVEPostProcess2D(self.job_name)

    def _create_part(self):
        PartGenerator = FullCircleInclusion(
            self.model,
            self.loc_info,
            self.model_name,
            self.part_name,
            self.instance_name,
        )
        PartGenerator.create_part()

    def _create_assembly(self):

        self.model.rootAssembly.features.changeKey(
            fromName="Final_Stuff-1", toName="Final_Stuff"
        )
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Final_Stuff"]

    def _create_geometry_set(self):

        # create part set for faces
        p = self.model.parts[self.part_name]
        faces = p.faces[:]
        p.Set(faces=faces, name="all_faces")
        fiberface = p.faces.getByBoundingCylinder(
            (self.loc[0][0], self.loc[0][1], 0.0),
            (self.loc[0][0], self.loc[0][1], 1.0),
            self.radius + 0.001 * self.radius,
        )
        p.Set(faces=fiberface, name="fiberface")
        for ii in range(1, len(self.loc)):
            fiberface_1 = p.faces.getByBoundingCylinder(
                (self.loc[ii][0], self.loc[ii][1], 0.0),
                (self.loc[ii][0], self.loc[ii][1], 1.0),
                self.radius + 0.001 * self.radius,
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

    def _creat_mesh(self, name_edges):

        import mesh

        elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(
            elemCode=CPE6M,
            elemLibrary=STANDARD,
            secondOrderAccuracy=ON,
            distortionControl=DEFAULT,
        )
        self.part.setElementType(
            regions=(self.part.faces[:],), elemTypes=(elemType1, elemType2)
        )
        self.part.seedPart(
            size=self.mesh_size, deviationFactor=0.2, minSizeFactor=0.9
        )
        for ii in range(len(name_edges)):
            self.part.seedEdgeBySize(
                edges=self.part.sets[name_edges[ii]].edges,
                size=self.mesh_size,
                deviationFactor=0.9,
                constraint=FIXED,
            )
        self.part.generateMesh()

    def _create_material(self, name_faces):
        # assign material property to matrix material
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Matrix",
            model=self.model,
            part=self.part,
            name_set=name_faces[0],
        )
        MaterialGenerator.CreateVonMisesPlasticMaterial(
            E=self.E_matrix,
            v=self.Pr_matrix,
            yield_criterion=self.yield_crtierion_matrix,
        )
        # assign material property to fiber material
        MaterialGenerator = AbaqusMaterialLib(
            name_mat="Fiber",
            model=self.model,
            part=self.part,
            name_set=name_faces[1],
        )
        MaterialGenerator.CreateElasticMaterial(
            E=self.E_fiber, v=self.Pr_fiber
        )

    def _create_step(
        self,
        initialInc=0.1,
        maxInc=1,
        maxNumInc=1000000,
        minInc=1e-20,
        timePeriod=20,
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
                "ETOTAL",
            ),
            timeInterval=0.01,
        )

    def _create_path_load(self):
        load_creater = Loading2D(self.model, self.assembly, self.instance_name)
        load_creater.create_path_load(
            self.loads, self.loads_path, self.time_period
        )
