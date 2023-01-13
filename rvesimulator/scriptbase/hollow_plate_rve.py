#                                                                       Modules
# =============================================================================
# standard
import sys

# third party
import numpy

# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *

try:
    import cPickle as pickle  # Improve speed
except ValueError:
    import pickle

# local
from base import RVE2DBase
from geometry import HollowPlate
from material import AbaqusMaterialLib
from postprocess import RVEPostProcess2D

#                                                          Authorship & Credits
# =============================================================================
__author__ = "Jiaxiang Yi (J.Yi@tudelft.nl)"
__credits__ = ["Jiaxiang Yi"]
__status__ = "Stable"
# =============================================================================
#
# =============================================================================


class HollowPlateRVE(RVE2DBase):
    def __init__(self, sim_info=None):
        """Initialization of the hollow plate rve simulation

        Parameters
        ----------
        sim_info : dict, optional
            simulation information, by default None
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
        """assemble the simulation"""

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
        if self.platform == "cluster":
            RVEPostProcess2D(self.job_name)

    def _create_part(self):
        """create part"""
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
        """create material

        Parameters
        ----------
        name_faces : str
            face name for that material
        """
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
        """create geometry set in abaqus

        Returns
        -------
        [list, list, list]
            sets of different geometry part
        """

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
        maxInc=1.0,
        maxNumInc=1000000,
        minInc=1e-20,
        timePeriod=1.0,
    ):
        """create step for the simulation

        Parameters
        ----------
        initialInc : float, optional
            initial increment, by default 0.1
        maxInc : float, optional
            maximum increment , by default 1.0
        maxNumInc : int, optional
            maximum number of incement, by default 1000000
        minInc :float, optional
             minimum incremnent step, by default 1e-20
        timePeriod : float, optional
            total simulation time , by default 1
        """

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
