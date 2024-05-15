# abaqus modulus
# system and standard python packages
import sys

import numpy
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *

try:
    import cPickle as pickle

except ValueError:
    import pickle

# import local functions
from basic_analysis_scripts.common_procedure import CommonProcedure
from basic_analysis_scripts.geometry_modeling import MultiCirclesPlates
from basic_analysis_scripts.jobs import Jobs
from basic_analysis_scripts.loading_condition import (
    HistoryDependentDisplacement2D, NormalDisplacementLoading)
from basic_analysis_scripts.materials import MaterialLib
from basic_analysis_scripts.mesh import Mesh2D
from basic_analysis_scripts.periodical_boundary_condition import \
    PeriodicalBoundaryCondition2D
from basic_analysis_scripts.post_process import PostProcess2D
from basic_analysis_scripts.steps import (LargeDeformationSteps,
                                          SmallDeformationSteps)


# =============================================================================
class SingleMaterialSVEBase:

    def create_pbc(self, vertices_name):
        # create reference points
        self.create_reference_points()
        # create pbc for vertices
        self.pbc_for_vertices(
            set_name=["NodeLB", "NodeLT", "NodeRB", "NodeRT"],
            geometry_name=vertices_name,
        )
        # create pbc for edges
        self.pbc_for_edges()

    def create_geometry_set(self):

        # create sets for faces:
        faces = self.model.parts[self.part_name].faces[:]
        self.part.Set(faces=faces, name="all_faces")
        faces_name = ["all_faces"]
        # create sets for edges:
        edges_name = self.create_sets_for_edges()
        # create sets for vertices
        vertices_name = self.create_sets_for_vertices()
        return faces_name, edges_name, vertices_name

# =============================================================================


class ElasticRegularLoads(
    CommonProcedure,
    MaterialLib,
    Mesh2D,
    MultiCirclesPlates,
    PeriodicalBoundaryCondition2D,
    NormalDisplacementLoading,
    SmallDeformationSteps,
    Jobs,
    SingleMaterialSVEBase
):

    def __init__(
        self,
        sim_info={
            "job_name": "single_material_sve",
            "location_information": None,
            "size": None,
            "youngs_modulus": None,
            "poisson_ratio": None,
            "mesh_partition": None,
            "strain": None,
            "num_cpu": None,
        },
    ):
        """
        Initialization function of Hollow Plate SVE case
        """
        # define the names of RVE simulation and set them as private variables
        self.model_name = "SVE"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = str(sim_info["job_name"])
        self.num_cpu = sim_info["num_cpu"]

        # define the import elements of RVE
        self.model = None
        self.sketch = None
        self.part = None
        self.assembly = None
        self.material = None

        #  define the geometry information
        self.length = sim_info["size"]
        self.width = sim_info["size"]
        self.circles_information = sim_info["location_information"]
        self.center = [self.length / 2.0, self.width / 2.0]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )

        # loading condition
        self.strain = sim_info["strain"]

        # material parameters
        self.youngs_modulus = sim_info["youngs_modulus"]
        self.poisson_ratio = sim_info["poisson_ratio"]

        # simulation time information
        self.time_period = 1.0
        self.time_interval = 0.1

        # execute the simulation
        self.script_generator()

    def script_generator(self):
        # create a new model
        self.create_new_model()
        # delete existed models (if any)
        self.delete_existed_models()
        # create part
        self.create_part()
        # create assembly for abaqus analysis
        self.create_assembly()
        # create geometry set for faces, edges, and vertices
        faces_name, edges_name, vertices_name = self.create_geometry_set()
        # create material
        self.create_elastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus,
            poisson=self.poisson_ratio,
        )
        # create mesh
        self.create_mesh(edges_name=edges_name, element_type="quadratic")
        # create pbc
        self.create_pbc(vertices_name=vertices_name)
        # create step
        self.create_step()
        # create loading
        self.create_load()
        # create job
        self.create_sequential_job(subroutine_path="")

        # # post process
        # if self.platform == "cluster" or self.platform == "windows":
        #     PostProcess2D(self.job_name)

# =============================================================================


class VonMisesPlasticRegularLoads(CommonProcedure,
                                  MaterialLib,
                                  Mesh2D,
                                  MultiCirclesPlates,
                                  PeriodicalBoundaryCondition2D,
                                  NormalDisplacementLoading,
                                  LargeDeformationSteps,
                                  Jobs,
                                  SingleMaterialSVEBase):

    def __init__(
        self,
        sim_info={
            "job_name": "single_material_sve",
            "location_information": None,
            "size": None,
            "youngs_modulus": None,
            "poisson_ratio": None,
            "hardening_table": None,
            "mesh_partition": None,
            "strain": None,
            "num_cpu": None,
            "num_steps": None,
            "simulation_time": None},
    ):
        """
        Initialization function of Hollow Plate SVE case
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
        self.length = sim_info["size"]
        self.width = sim_info["size"]
        self.center = [self.length / 2.0, self.width / 2.0]
        self.circles_information = sim_info["location_information"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )

        # loading condition
        self.strain = sim_info["strain"]

        # simulation time information
        self.time_period = sim_info["simulation_time"]
        self.time_interval = sim_info["simulation_time"] / \
            sim_info["num_steps"]

        # material properties
        self.youngs_modulus = sim_info["youngs_modulus"]
        self.poisson_ratio = sim_info["poisson_ratio"]
        self.hardening_table = numpy.zeros(
            (len(sim_info["hardening_table"][0]), 2)
        )
        for ii in range(2):
            self.hardening_table[:, ii] = sim_info["hardening_table"][ii]
        # execute the simulation
        self.script_generator()

    def script_generator(self):
        # create a new model
        self.create_new_model()
        # delete existed models (if any)
        self.delete_existed_models()
        # create part
        self.create_part()
        # create assembly for abaqus analysis
        self.create_assembly()
        # create geometry set for faces, edges, and vertices
        faces_name, edges_name, vertices_name = self.create_geometry_set()
        # create material
        self.create_von_mises_plastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus,
            poisson=self.poisson_ratio,
            hardening_table=self.hardening_table,
        )
        # create mesh
        self.create_mesh(edges_name=edges_name, element_type="quadratic")
        # create pbc
        self.create_pbc(vertices_name=vertices_name)
        # create step
        self.create_step()
        # create loading
        self.create_load()
        # create job
        self.create_sequential_job(subroutine_path="")

        # # post process
        # if self.platform == "cluster" or self.platform == "windows":
        #     PostProcess2D(self.job_name)

# =============================================================================


class VonMisesPlasticPathLoads(CommonProcedure,
                               MaterialLib,
                               Mesh2D,
                               MultiCirclesPlates,
                               PeriodicalBoundaryCondition2D,
                               HistoryDependentDisplacement2D,
                               LargeDeformationSteps,
                               Jobs,
                               SingleMaterialSVEBase):
    def __init__(
        self,
        sim_info={
            "job_name": "single_material_sve",
            "location_information": None,
            "size": None,
            "youngs_modulus": None,
            "poisson_ratio": None,
            "hardening_table": None,
            "mesh_partition": None,
            "strain": None,
            "strain_amplitude": None,
            "num_cpu": None,
            "num_steps": None,
            "simulation_time": None},
    ):
        """
        Initialization function of Hollow Plate SVE case
        """
        # define the names of RVE simulation and set them as private variables
        self.model_name = "SVE"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = str(sim_info["job_name"])
        self.num_cpu = sim_info["num_cpu"]

        # define the import elements of RVE
        self.model = None
        self.sketch = None
        self.part = None
        self.assembly = None
        self.material = None

        #  define the geometry information
        self.length = sim_info["size"]
        self.width = sim_info["size"]
        self.center = [self.length / 2.0, self.width / 2.0]
        self.circles_information = sim_info["location_information"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )

        # loading condition
        # strain
        self.strain = sim_info["strain"]

        # strain amplitude
        self.strain_amplitude = numpy.zeros(
            (len(sim_info["strain_amplitude"][0]), 3)
        )
        for ii in range(3):
            self.strain_amplitude[:, ii] = sim_info["strain_amplitude"][ii]

        # simulation time information
        self.time_period = sim_info["simulation_time"]
        self.time_interval = sim_info["simulation_time"] / \
            sim_info["num_steps"]

        # material properties
        self.youngs_modulus = sim_info["youngs_modulus"]
        self.poisson_ratio = sim_info["poisson_ratio"]
        self.hardening_table = numpy.zeros(
            (len(sim_info["hardening_table"][0]), 2)
        )
        for ii in range(2):
            self.hardening_table[:, ii] = sim_info["hardening_table"][ii]

        # execute the simulation
        self.script_generator()

    def script_generator(self):
        # create a new model
        self.create_new_model()
        # delete existed models (if any)
        self.delete_existed_models()
        # create part
        self.create_part()
        # create assembly for abaqus analysis
        self.create_assembly()
        # create geometry set for faces, edges, and vertices
        faces_name, edges_name, vertices_name = self.create_geometry_set()
        # create material
        self.create_von_mises_plastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus,
            poisson=self.poisson_ratio,
            hardening_table=self.hardening_table,
        )
        # create mesh
        self.create_mesh(edges_name=edges_name, element_type="quadratic")
        # create pbc
        self.create_pbc(vertices_name=vertices_name)
        # create step
        self.create_step()
        # create loading
        self.create_load()
        # create job
        self.create_sequential_job(subroutine_path="")

        # # post process
        # if self.platform == "cluster" or self.platform == "windows":
        #     PostProcess2D(self.job_name)
