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

except:
    import pickle

# import local functions
from basic_analysis_scripts.common_procedure import CommonProcedure
from basic_analysis_scripts.geometry_modeling import MultiCirclesInclusion
from basic_analysis_scripts.jobs import Jobs
from basic_analysis_scripts.loading_condition import (
    HistoryDependentDisplacement2D, NormalDisplacementLoading)
from basic_analysis_scripts.materials import MaterialLib
from basic_analysis_scripts.mesh import Mesh2D
from basic_analysis_scripts.periodical_boundary_condition import \
    PeriodicalBoundaryCondition2D
from basic_analysis_scripts.post_process import PostProcess2D
from basic_analysis_scripts.steps import LargeDeformationSteps


# =============================================================================
class TwoMaterialsRVEBase(CommonProcedure,
                          MaterialLib,
                          Mesh2D,
                          PeriodicalBoundaryCondition2D,
                          Jobs,
                          LargeDeformationSteps,
                          MultiCirclesInclusion):

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

    def create_assembly(self):

        self.model.rootAssembly.features.changeKey(
            fromName="Final_Stuff-1", toName="Final_Stuff"
        )
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Final_Stuff"]

    def create_geometry_set(self):

        # get faces
        faces = self.model.parts[self.part_name].faces[:]
        # for all faces
        self.part.Set(faces=faces, name="all_faces")

        fiberface = self.part.faces.getByBoundingBox(
            xMin=self.circles_information[0][0] -
            self.circles_information[0][2] - 0.001 *
            self.circles_information[0][2],
            xMax=self.circles_information[0][0] +
            self.circles_information[0][2] + 0.001 *
            self.circles_information[0][2],
            yMin=self.circles_information[0][1] -
            self.circles_information[0][2] - 0.001 *
            self.circles_information[0][2],
            yMax=self.circles_information[0][1] +
            self.circles_information[0][2] + 0.001 *
            self.circles_information[0][2],
            zMin=0.0,
            zMax=1.0,
        )

        # fiber faces
        self.part.Set(faces=fiberface, name="fiberface")
        for ii in range(1, len(self.circles_information)):

            fiberface_1 = self.part.faces.getByBoundingBox(
                xMin=self.circles_information[ii][0] -
                self.circles_information[ii][2] - 0.001 *
                self.circles_information[ii][2],
                xMax=self.circles_information[ii][0] +
                self.circles_information[ii][2] + 0.001 *
                self.circles_information[ii][2],
                yMin=self.circles_information[ii][1] -
                self.circles_information[ii][2] - 0.001 *
                self.circles_information[ii][2],
                yMax=self.circles_information[ii][1] +
                self.circles_information[ii][2] + 0.001 *
                self.circles_information[ii][2],
                zMin=0.0,
                zMax=1.0,
            )
            self.part.Set(faces=fiberface_1, name="fiberface_1")
            self.part.SetByBoolean(
                name="fiberface",
                sets=(self.part.sets["fiberface_1"],
                      self.part.sets["fiberface"]),
                operation=UNION,
            )

        # delete fiberface 1
        del self.part.sets["fiberface_1"]
        self.part.SetByBoolean(
            name="matrixface",
            sets=(self.part.sets["all_faces"], self.part.sets["fiberface"]),
            operation=DIFFERENCE,
        )
        faces_name = ["matrixface", "fiberface"]

        # create sets for edges:
        edges_name = self.create_sets_for_edges()
        vertices_name = self.create_sets_for_vertices()

        return faces_name, edges_name, vertices_name


# =============================================================================
class VonMisesPlasticElasticRegularLoads(NormalDisplacementLoading,
                                         TwoMaterialsRVEBase):

    def __init__(
        self,
        sim_info={
            "job_name": "two_materials_rve",
            "location_information": None,
            "len_start": None,
            "len_end": None,
            "wid_start": None,
            "wid_end": None,
            "radius_mu": None,
            "radius_std": None,
            "youngs_modulus_matrix": None,
            "poisson_ratio_matrix": None,
            "hardening_table": None,
            "youngs_modulus_fiber": None,
            "poisson_ratio_fiber": None,
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
        self.model_name = "RVE"
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
        self.len_start = sim_info["len_start"]
        self.len_end = sim_info["len_end"]
        self.wid_start = sim_info["wid_start"]
        self.wid_end = sim_info["wid_end"]
        self.circles_information = sim_info["location_information"]
        # modeling
        self.length = (sim_info["len_end"] - sim_info["len_start"]
                       ) - 2 * sim_info["radius_mu"]

        self.width = (sim_info["wid_end"] - sim_info["wid_start"]
                      ) - 2 * sim_info["radius_mu"]

        self.center = [
            (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
            (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
        ]
        self.radius_mu = sim_info["radius_mu"]
        self.radius_std = sim_info["radius_std"]
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
        self.youngs_modulus_fiber = sim_info["youngs_modulus_fiber"]
        self.poisson_ratio_fiber = sim_info["poisson_ratio_fiber"]
        self.youngs_modulus_matrix = sim_info["youngs_modulus_matrix"]
        self.poisson_ratio_matrix = sim_info["poisson_ratio_matrix"]
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
        # create material for matrix material
        self.create_von_mises_plastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus_matrix,
            poisson=self.poisson_ratio_matrix,
            hardening_table=self.hardening_table,
        )
        # create material for fiber material
        self.create_elastic_material(
            geometry_set_name=faces_name[1],
            material_name="fiber",
            youngs=self.youngs_modulus_fiber,
            poisson=self.poisson_ratio_fiber,

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

class TwoElasticRegularLoads(NormalDisplacementLoading,
                            TwoMaterialsRVEBase):

    def __init__(
        self,
        sim_info={
            "job_name": "elastic_rve",
            "location_information": None,
            "len_start": None,
            "len_end": None,
            "wid_start": None,
            "wid_end": None,
            "radius_mu": None,
            "radius_std": None,
            "youngs_modulus_matrix": None,
            "poisson_ratio_matrix": None,
            "youngs_modulus_fiber": None,
            "poisson_ratio_fiber": None,
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
        self.model_name = "RVE"
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
        self.len_start = sim_info["len_start"]
        self.len_end = sim_info["len_end"]
        self.wid_start = sim_info["wid_start"]
        self.wid_end = sim_info["wid_end"]
        self.circles_information = sim_info["location_information"]
        # modeling
        self.length = (sim_info["len_end"] - sim_info["len_start"]
                       ) - 2 * sim_info["radius_mu"]

        self.width = (sim_info["wid_end"] - sim_info["wid_start"]
                      ) - 2 * sim_info["radius_mu"]

        self.center = [
            (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
            (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
        ]
        self.radius_mu = sim_info["radius_mu"]
        self.radius_std = sim_info["radius_std"]
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
        self.youngs_modulus_fiber = sim_info["youngs_modulus_fiber"]
        self.poisson_ratio_fiber = sim_info["poisson_ratio_fiber"]
        self.youngs_modulus_matrix = sim_info["youngs_modulus_matrix"]
        self.poisson_ratio_matrix = sim_info["poisson_ratio_matrix"]

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
        # create material for matrix material
        self.create_elastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus_matrix,
            poisson=self.poisson_ratio_matrix,
        )
        # create material for fiber material
        self.create_elastic_material(
            geometry_set_name=faces_name[1],
            material_name="fiber",
            youngs=self.youngs_modulus_fiber,
            poisson=self.poisson_ratio_fiber,

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


# =============================================================================


class VonMisesPlasticElasticPathLoads(HistoryDependentDisplacement2D,
                                      TwoMaterialsRVEBase):
    def __init__(
        self,
        sim_info={
            "job_name": "single_material_rve",
            "location_information": None,
            "len_start": None,
            "len_end": None,
            "wid_start": None,
            "wid_end": None,
            "radius_mu": None,
            "radius_std": None,
            "youngs_modulus_matrix": None,
            "poisson_ratio_matrix": None,
            "hardening_table": None,
            "youngs_modulus_fiber": None,
            "poisson_ratio_fiber": None,
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
        self.model_name = "RVE"
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
        self.len_start = sim_info["len_start"]
        self.len_end = sim_info["len_end"]
        self.wid_start = sim_info["wid_start"]
        self.wid_end = sim_info["wid_end"]
        self.circles_information = sim_info["location_information"]
        # modeling
        self.length = (sim_info["len_end"] - sim_info["len_start"]
                       ) - 2 * sim_info["radius_mu"]

        self.width = (sim_info["wid_end"] - sim_info["wid_start"]
                      ) - 2 * sim_info["radius_mu"]

        self.center = [
            (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
            (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
        ]
        self.radius_mu = sim_info["radius_mu"]
        self.radius_std = sim_info["radius_std"]
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
        self.youngs_modulus_fiber = sim_info["youngs_modulus_fiber"]
        self.poisson_ratio_fiber = sim_info["poisson_ratio_fiber"]
        self.youngs_modulus_matrix = sim_info["youngs_modulus_matrix"]
        self.poisson_ratio_matrix = sim_info["poisson_ratio_matrix"]
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
        # create material for matrix material
        self.create_von_mises_plastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus_matrix,
            poisson=self.poisson_ratio_matrix,
            hardening_table=self.hardening_table,
        )
        # create material for fiber material
        self.create_elastic_material(
            geometry_set_name=faces_name[1],
            material_name="fiber",
            youngs=self.youngs_modulus_fiber,
            poisson=self.poisson_ratio_fiber,

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



# =============================================================================
class VonMisesPlasticMatrixFiberRegularLoads(NormalDisplacementLoading,
                                             TwoMaterialsRVEBase):

    def __init__(
        self,
        sim_info={
            "job_name": "plastic_rve",
            "location_information": None,
            "len_start": None,
            "len_end": None,
            "wid_start": None,
            "wid_end": None,
            "radius_mu": None,
            "radius_std": None,
            "youngs_modulus_matrix": None,
            "poisson_ratio_matrix": None,
            "hardening_table_matrix": None,
            "youngs_modulus_fiber": None,
            "poisson_ratio_fiber": None,
            "hardening_table_fiber": None,
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
        self.model_name = "RVE"
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
        self.len_start = sim_info["len_start"]
        self.len_end = sim_info["len_end"]
        self.wid_start = sim_info["wid_start"]
        self.wid_end = sim_info["wid_end"]
        self.circles_information = sim_info["location_information"]
        # modeling
        self.length = (sim_info["len_end"] - sim_info["len_start"]
                       ) - 2 * sim_info["radius_mu"]

        self.width = (sim_info["wid_end"] - sim_info["wid_start"]
                      ) - 2 * sim_info["radius_mu"]

        self.center = [
            (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
            (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
        ]
        self.radius_mu = sim_info["radius_mu"]
        self.radius_std = sim_info["radius_std"]
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
        self.youngs_modulus_fiber = sim_info["youngs_modulus_fiber"]
        self.poisson_ratio_fiber = sim_info["poisson_ratio_fiber"]
        self.hardening_table_fiber = numpy.zeros(
            (len(sim_info["hardening_table_fiber"][0]), 2)
        )
        for ii in range(2):
            self.hardening_table_fiber[:, ii] = sim_info["hardening_table_fiber"][ii]

        self.youngs_modulus_matrix = sim_info["youngs_modulus_matrix"]
        self.poisson_ratio_matrix = sim_info["poisson_ratio_matrix"]
        self.hardening_table_matrix = numpy.zeros(
            (len(sim_info["hardening_table_matrix"][0]), 2))
        for ii in range(2):
            self.hardening_table_matrix[:, ii] = sim_info["hardening_table_matrix"][ii]

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
        # create material for matrix material
        self.create_von_mises_plastic_material(
            geometry_set_name=faces_name[0],
            material_name="matrix",
            youngs=self.youngs_modulus_matrix,
            poisson=self.poisson_ratio_matrix,
            hardening_table=self.hardening_table_matrix,
        )
        # create material for fiber material
        self.create_von_mises_plastic_material(
            geometry_set_name=faces_name[1],
            material_name="fiber",
            youngs=self.youngs_modulus_fiber,
            poisson=self.poisson_ratio_fiber,
            hardening_table=self.hardening_table_fiber,
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
