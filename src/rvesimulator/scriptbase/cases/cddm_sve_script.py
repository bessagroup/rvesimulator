# extend the system path
import sys

# import python libraries
import numpy

# abaqus
from abaqus import *
from abaqusConstants import *

# local function
from basic_script.base import RVE2DBase
from basic_script.geometry import HeterHollowPlate
from basic_script.loadings import Loading2D
from basic_script.material import AbaqusMaterialLib
from basic_script.postprocess import RVEPostProcess2D
from caeModules import *


class CDDMSVE(RVE2DBase):
    def __init__(
        self,
        sim_info={
            "length": 1.0,
            "width": 1.0,
            "location_information": None,
            "youngs_modulus": None,
            "poisson_ratio": None,
            "hardening_table": None,
            "mesh_partition": None,
            "strain": [0.02, 0.02, 0.02],
            "strain_amplitude": None,
            "simulation_time": 1.0,
            "num_steps": 100,
            "job_name": "cddm_sve",
            "num_cpu": 1,
            "platform": "ubuntu",
        },
    ):
        """
        Initializaton function of CDDM_SVE case

        """

        # define the names of RVE simulation and set them as private variables
        self.model_name = "SVE"
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
        self.strain = sim_info["strain"]

        # material properties
        self.youngs_modulus = sim_info["youngs_modulus"]
        self.poisson_ratio = sim_info["poisson_ratio"]
        # simulation time and time intervel for recording the results
        self.time_period = sim_info["simulation_time"]
        self.time_interval = (
            sim_info["simulation_time"] / sim_info["num_steps"]
        )
        # get the hardening law
        self.hardening_table = numpy.zeros(
            (len(sim_info["hardening_table"][0]), 2)
        )
        for ii in range(2):
            self.hardening_table[:, ii] = sim_info["hardening_table"][ii]

        # get the strain amplitude
        self.strain_amplitude = numpy.zeros(
            (len(sim_info["strain_amplitude"][0]), 3)
        )
        for ii in range(3):
            self.strain_amplitude[:, ii] = sim_info["strain_amplitude"][ii]

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
        self._create_step()
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
            v=self.poisson_ratio,
            yield_criterion=self.hardening_table,
        )

    def _create_path_load(self):
        load_creater = Loading2D(self.model, self.assembly, self.instance_name)
        load_creater.create_path_load(
            self.strain, self.strain_amplitude, self.time_period
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
    ):

        self.model.StaticStep(name="Step-1", previous="Initial", nlgeom=ON)
        self.model.StaticStep(
            initialInc=initialInc,
            maxInc=maxInc,
            maxNumInc=maxNumInc,
            minInc=minInc,
            name="Step-1",
            previous="Initial",
            timePeriod=self.time_period,
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
            ),
            timeInterval=self.time_interval,
        )
