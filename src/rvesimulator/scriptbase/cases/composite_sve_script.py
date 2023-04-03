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

from scriptbase.basic_script.base import RVE2DBase
from scriptbase.basic_script.geometry import FullCircleInclusion
from scriptbase.basic_script.loadings import Loading2D
from scriptbase.basic_script.material import AbaqusMaterialLib
from scriptbase.basic_script.postprocess import RVEPostProcess2D


class SVE(RVE2DBase):
    def __init__(
        self,
        sim_info={
            "location_information": None,
            "length": None,
            "width": None,
            "job_name": "sve",
            "strain": [0.05, 0.0, 0.0],
            "strain_amplitude": None,
            "E_matrix": None,
            "Pr_matrix": None,
            "hardening_table": None,
            "E_fiber": None,
            "Pr_fiber": None,
            "mesh_partition": 50,
            "simulation_time": 1.0,
            "num_steps": None,
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
        self.length = sim_info["length"]
        self.width = sim_info["width"]
        self.center = [self.length / 2.0, self.width / 2.0]

        # information of RVE modeling
        self.strain = sim_info["strain"]
        # self.loads_path = sim_info["loads_path"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )
        # simulation time and time intervel for recording the results
        self.time_period = sim_info["simulation_time"]
        self.time_interval = (
            sim_info["simulation_time"] / sim_info["num_steps"]
        )
        # material properties
        self.E_matrix = sim_info["E_matrix"]
        self.Pr_matrix = sim_info["Pr_matrix"]
        self.E_fiber = sim_info["E_fiber"]
        self.Pr_fiber = sim_info["Pr_fiber"]

        # opreation on transfer list to numpy
        self.hardening_table = numpy.zeros(
            (len(sim_info["hardening_table"][0]), 2)
        )
        for ii in range(2):
            self.hardening_table[:, ii] = sim_info["hardening_table"][ii]

        self.strain_amplitude = numpy.zeros(
            (len(sim_info["strain_amplitude"][0]), 3)
        )
        for ii in range(3):
            self.strain_amplitude[:, ii] = sim_info["strain_amplitude"][ii]

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
        self._create_step()
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
            yield_criterion=self.hardening_table,
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
    ):
        self.model.StaticStep(name="Step-1", previous="Initial")
        self.model.StaticStep(
            initialInc=initialInc,
            maxInc=maxInc,
            maxNumInc=maxNumInc,
            minInc=minInc,
            name="Step-1",
            previous="Initial",
            timePeriod=self.time_period,
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
            timeInterval=self.time_interval,
        )
        self.model.FieldOutputRequest(
            name="F-Output-2",
            createStepName="Step-1",
            variables=("U", "RF"),
            timeInterval=self.time_interval,
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
            timeInterval=self.time_interval,
        )

    def _create_path_load(self):
        load_creater = Loading2D(self.model, self.assembly, self.instance_name)
        load_creater.create_path_load(
            self.strain, self.strain_amplitude, self.time_period
        )
