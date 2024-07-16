
# -*- coding: mbcs -*-
import assembly
import mesh
import numpy
import regionToolset
from math import ceil
from abaqus import *
from abaqusConstants import *
from caeModules import *
# import packages for abaqus post-processing
from odbAccess import *

try:
    import cPickle as pickle  # Improve speed
except:
    import pickle


class PPPEMixtureCohesive:
    """ PP and PE 2D RVE with cohesive elements"""

    def __init__(self, sim_info={"location_information": None,
                                 "radius_mu": None,
                                 "radius_std": None,
                                 "len_start": None,
                                 "len_end": None,
                                 "wid_start": None,
                                 "wid_end": None,
                                 "job_name": "pp_pe_cohesive_2d",
                                 "strain": [0.10, 0.0, 0.0],
                                 "paras_pp": None,
                                 "paras_pe": None,
                                 "damage_stress": None,
                                 "damage_energy": None,
                                 "power_law_exponent_cohesive": 1.0,
                                 "mesh_partition": None,
                                 "simulation_time": 100.0,
                                 "num_steps": None,
                                 "record_time_step": 100,
                                 "young_modulus_cohesive": None,
                                 "num_cpu": 1,
                                 "radius_cohesive_factor": 1.02,
                                 "damage_onset_criteria": "MaxStress",
                                 "subroutine_path": None}):
        # ------------------------------ parameters ------------------------- #
        # names of model, part, instance
        self.model_name = "pp_pe_cohesive_2d"
        self.part_name = "Final_Stuff"
        self.instance_name = "Final_Stuff"
        self.job_name = str(sim_info["job_name"])
        self.num_cpu = sim_info["num_cpu"]
        self.record_time_step = sim_info["record_time_step"]
        # information of geometry of RVE
        self.loc_info = sim_info
        # mech and sets information
        self.circles_information = sim_info["location_information"]
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
        self.strain = sim_info["strain"]
        self.mesh_size = (
            min(self.length, self.width) / sim_info["mesh_partition"]
        )
        self.time_period = sim_info["simulation_time"]
        self.time_interval = (
            sim_info["simulation_time"] / sim_info["num_steps"]
        )

        # material properties
        self.paras_pp = sim_info["paras_pp"]
        self.paras_pe = sim_info["paras_pe"]
        self.subroutine_path = str(sim_info["subroutine_path"])

        # damage parameters
        self.damage_onset_criteria = sim_info["damage_onset_criteria"]
        self.damage_stress = sim_info["damage_stress"]
        self.damage_energy = sim_info["damage_energy"]
        self.young_modulus_cohesive = sim_info["young_modulus_cohesive"]
        self.cohesive_radius_factor = sim_info["radius_cohesive_factor"]
        self.power = sim_info["power_law_exponent_cohesive"]
        # submit
        self.create_simulation_job()
        self.create_job()
        # self.submit_job()

        # # post process
        # if sim_info["platform"] == "cluster" or sim_info["platform"] == "windows":
        #     # post process for getting the results
        #     PostProcess(job_name=self.job_name,
        #                 record_time_step=self.record_time_step)

    def create_simulation_job(self):

        # get delta
        delta = min(min(self.length, self.width) / 1000, self.mesh_size / 10)

        # get model information
        instance_name = self.instance_name
        model_name = self.model_name

        # create model---------------------------------------------------------
        mdb.models.changeKey(fromName='Model-1', toName=model_name)
        model = mdb.models[model_name]
        # delete existed model
        if 'Model-1' in mdb.models.keys():
            del mdb.models['Model-1']

        # create sketch--------------------------------------------------------
        # model the matrix part first
        sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        sketch.rectangle(point1=(0.0, 0.0), point2=(self.length, self.width))
        part = model.Part(name='Final_Stuff',
                          dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
        part.BaseShell(sketch=sketch)
        del sketch

        # partition of square
        faces = part.faces
        t = part.MakeSketchTransform(
            sketchPlane=faces[0], sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
        # origin if used for shifting the coordinates
        sketch = model.ConstrainedSketch(
            name='__profile__', sheetSize=2.82, gridSpacing=0.07, transform=t)
        sketch.setPrimaryObject(option=SUPERIMPOSE)
        part.projectReferencesOntoSketch(sketch=sketch, filter=COPLANAR_EDGES)

        # create the fibers
        for ii in range(len(self.circles_information)):
            # actual fibers
            sketch.CircleByCenterPerimeter(
                center=(self.circles_information[ii][0],
                        self.circles_information[ii][1]),
                point1=(self.circles_information[ii][0] +
                        self.circles_information[ii][2], self.circles_information[ii][1]),
            )
            # cohesive zone
            sketch.CircleByCenterPerimeter(
                center=(self.circles_information[ii][0],
                        self.circles_information[ii][1]),
                point1=(self.circles_information[ii][0] + self.cohesive_radius_factor *
                        self.circles_information[ii][2], self.circles_information[ii][1]),
            )
        part.PartitionFaceBySketch(faces=faces[:], sketch=sketch)
        del sketch

        # create sets----------------------------------------------------------
        # get faces
        faces = part.faces[:]

        # for all faces
        part.Set(faces=faces, name="all_faces")

        # fiber faces
        fiberface = part.faces.getByBoundingCylinder(
            (self.circles_information[0][0],
             self.circles_information[0][1], 0.0),
            (self.circles_information[0][0],
             self.circles_information[0][1], 1.0),
            self.circles_information[0][2] + 0.001 *
            self.circles_information[0][2],
        )
        part.Set(faces=fiberface, name="fiberface")
        for ii in range(1, len(self.circles_information)):
            fiberface_1 = part.faces.getByBoundingCylinder(
                (self.circles_information[ii][0],
                 self.circles_information[ii][1], 0.0),
                (self.circles_information[ii][0],
                 self.circles_information[ii][1], 1.0),
                self.circles_information[ii][2] +
                0.001 * self.circles_information[ii][2],
            )
            part.Set(faces=fiberface_1, name="fiberface_1")
            part.SetByBoolean(
                name="fiberface",
                sets=(part.sets["fiberface_1"],  part.sets["fiberface"]),
                operation=UNION,
            )
        # delete fiber face 1
        del part.sets["fiberface_1"]

        # create set for cohesive and fibers
        cohesive_fiber = part.faces.getByBoundingCylinder(
            (self.circles_information[0][0],
             self.circles_information[0][1], 0.0),
            (self.circles_information[0][0],
             self.circles_information[0][1], 1.0),
            self.circles_information[0][2]*self.cohesive_radius_factor +
            0.001 * self.circles_information[0][2],
        )

        part.Set(faces=cohesive_fiber, name="cohesive_fiber")
        for ii in range(1, len(self.circles_information)):
            cohesive_fiber_1 = part.faces.getByBoundingCylinder(
                (self.circles_information[ii][0],
                 self.circles_information[ii][1], 0.0),
                (self.circles_information[ii][0],
                 self.circles_information[ii][1], 1.0),
                self.circles_information[ii][2]*self.cohesive_radius_factor +
                0.001 * self.circles_information[ii][2],
            )
            part.Set(faces=cohesive_fiber_1, name="cohesive_fiber_1")
            part.SetByBoolean(
                name="cohesive_fiber",
                sets=(part.sets["cohesive_fiber_1"],
                      part.sets["cohesive_fiber"]),
                operation=UNION,
            )
        # delete fiber face 1
        del part.sets["cohesive_fiber_1"]

        # create set for matrix
        part.SetByBoolean(
            name="matrixface",
            sets=(part.sets["all_faces"], part.sets["cohesive_fiber"]),
            operation=DIFFERENCE,
        )

        part.SetByBoolean(
            name="cohesive_face",
            sets=(part.sets["cohesive_fiber"], part.sets["fiberface"]),
            operation=DIFFERENCE,
        )

        del part.sets["cohesive_fiber"]

        # create sets for edges
        s = part.edges
        edgesLEFT = s.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] - self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        part.Set(edges=edgesLEFT, name="edgesLEFT")
        edgesRIGHT = s.getByBoundingBox(
            self.center[0] + self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        part.Set(edges=edgesRIGHT, name="edgesRIGHT")
        edgesTOP = s.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] + self.width / 2 - delta,
            0,
            self.center[0] + self.width / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        part.Set(edges=edgesTOP, name="edgesTOP")
        edgesBOT = s.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] - self.width / 2 + delta,
            0,
        )
        part.Set(edges=edgesBOT, name="edgesBOT")
        name_edges = ["edgesLEFT", "edgesRIGHT", "edgesTOP", "edgesBOT"]

        # create set for cohesive edges
        cohesive_edges = s.getByBoundingCylinder(
            (self.circles_information[0][0],
             self.circles_information[0][1], 0.0),
            (self.circles_information[0][0],
             self.circles_information[0][1], 1.0),
            self.circles_information[0][2]*self.cohesive_radius_factor +
            0.001 * self.circles_information[0][2],
        )
        part.Set(edges=cohesive_edges, name="cohesive_edges")

        for ii in range(1, len(self.circles_information)):
            cohesive_edges_1 = s.getByBoundingCylinder(
                (self.circles_information[ii][0],
                 self.circles_information[ii][1], 0.0),
                (self.circles_information[ii][0],
                 self.circles_information[ii][1], 1.0),
                self.circles_information[ii][2]*self.cohesive_radius_factor +
                0.001 * self.circles_information[ii][2],
            )
            part.Set(edges=cohesive_edges_1, name="cohesive_edges_1")
            part.SetByBoolean(
                name="cohesive_edges",
                sets=(part.sets["cohesive_edges_1"],
                      part.sets["cohesive_edges"]),
                operation=UNION,
            )
        # delete fiber face 1
        del part.sets["cohesive_edges_1"]

        # boolean opration to delete edges belong to matrix
        for ii in range(len(name_edges)):
            part.SetByBoolean(
                name="cohesive_edges",
                sets=(part.sets["cohesive_edges"],  part.sets[name_edges[ii]]),
                operation=DIFFERENCE,
            )

        # create set for vertices
        v = part.vertices
        vertexLB = v.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] - self.length / 2 + delta,
            self.center[1] - self.width / 2 + delta,
            0,
        )
        part.Set(vertices=vertexLB, name="VertexLB")
        vertexRB = v.getByBoundingBox(
            self.center[0] + self.length / 2 - delta,
            self.center[1] - self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] - self.width / 2 + delta,
            0,
        )
        part.Set(vertices=vertexRB, name="VertexRB")
        vertexRT = v.getByBoundingBox(
            self.center[0] + self.length / 2 - delta,
            self.center[1] + self.width / 2 - delta,
            0,
            self.center[0] + self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        part.Set(vertices=vertexRT, name="VertexRT")
        vertexLT = v.getByBoundingBox(
            self.center[0] - self.length / 2 - delta,
            self.center[1] + self.width / 2 - delta,
            0,
            self.center[0] - self.length / 2 + delta,
            self.center[1] + self.width / 2 + delta,
            0,
        )
        part.Set(vertices=vertexLT, name="VertexLT")
        vertices_name = ["VertexLB", "VertexLT", "VertexRB", "VertexRT"]

        # meshing ---------------------------------------------------------------------
        # calculate

        part.seedPart(deviationFactor=0.1,
                      minSizeFactor=0.05, size=self.mesh_size)
        # part.seedEdgeBySize(constraint=FINER, edges=part.sets['cohesive_edges'].edges,size=self.mesh_size/2.0)

        # mesh control
        part.setMeshControls(
            elemShape=QUAD, regions=part.sets['cohesive_face'].faces, technique=SWEEP)
        #

        elemType1 = mesh.ElemType(
            elemCode=COH2D4, elemLibrary=STANDARD, viscosity=0.0001)
        elemType2 = mesh.ElemType(elemCode=UNKNOWN_TRI, elemLibrary=STANDARD)

        part.setElementType(
            regions=part.sets["cohesive_face"], elemTypes=(
                elemType1, elemType2)
        )

        # element type (plane strain element)
        elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD,
                                  secondOrderAccuracy=OFF, hourglassControl=ENHANCED,
                                  distortionControl=DEFAULT)
        elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)
        part.setElementType(
            regions=part.sets["fiberface"], elemTypes=(elemType1, elemType2))
        part.setElementType(
            regions=part.sets["matrixface"], elemTypes=(elemType1, elemType2))

        part.generateMesh()

        # material properties ---------------------------------------------------------
        # material property for top part
        material_fiber = model.Material(name="material_fiber")
        material_fiber.Depvar(n=44)
        material_fiber.Density(table=((7.9e-9, ), ))
        material_fiber.UserMaterial(
            mechanicalConstants=(
                self.paras_pp[0],
                self.paras_pp[1],
                self.paras_pp[2],
                self.paras_pp[3],
                self.paras_pp[4],
                self.paras_pp[5],
                self.paras_pp[6],
                self.paras_pp[7],
                self.paras_pp[8],
                self.paras_pp[9],
                self.paras_pp[10],
                self.paras_pp[11],
                self.paras_pp[12],
                self.paras_pp[13],
            )
        )
        material_matrix = model.Material(name="material_matrix")
        material_matrix.Depvar(n=44)
        material_matrix.Density(table=((7.9e-9, ), ))
        material_matrix.UserMaterial(
            mechanicalConstants=(
                self.paras_pe[0],
                self.paras_pe[1],
                self.paras_pe[2],
                self.paras_pe[3],
                self.paras_pe[4],
                self.paras_pe[5],
                self.paras_pe[6],
                self.paras_pe[7],
                self.paras_pe[8],
                self.paras_pe[9],
                self.paras_pe[10],
                self.paras_pe[11],
                self.paras_pe[12],
                self.paras_pe[13],
            )
        )
        # cohesive material (parameters needs to be tuned)
        material_cohesive = model.Material(name='material_cohesive')
        material_cohesive.Density(table=((7.9e-9, ), ))
        material_cohesive.Elastic(type=TRACTION, table=((
            self.young_modulus_cohesive, self.young_modulus_cohesive, 0.0), ))
        if self.damage_onset_criteria == "MaxStress":
            material_cohesive.MaxsDamageInitiation(table=((
                self.damage_stress, self.damage_stress, 0.0), ))
            material_cohesive.maxsDamageInitiation.DamageEvolution(
                type=ENERGY, mixedModeBehavior=POWER_LAW,
                power=self.power,
                table=((self.damage_energy, self.damage_energy, 0.0), ))
            material_cohesive.maxsDamageInitiation.DamageStabilizationCohesive(
                cohesiveCoeff=0.000001)

        elif self.damage_onset_criteria == "QuadStress":
            material_cohesive.QuadsDamageInitiation(table=((
                self.damage_stress, self.damage_stress, 0.0), ))

            material_cohesive.quadsDamageInitiation.DamageEvolution(
                type=ENERGY, mixedModeBehavior=POWER_LAW,
                power=self.power,
                table=((self.damage_energy, self.damage_energy, 0.0), ))
            material_cohesive.quadsDamageInitiation.DamageStabilizationCohesive(
                cohesiveCoeff=0.000001)

        # create section and assign material property to corresponding section
        # matrix material
        model.HomogeneousSolidSection(
            name='matrix', material='material_matrix', thickness=None)
        part.SectionAssignment(region=part.sets['matrixface'], sectionName='matrix', offset=0.0,
                               offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

        # fiber material
        model.HomogeneousSolidSection(
            name='fiber', material='material_fiber', thickness=None)
        part.SectionAssignment(region=part.sets['fiberface'], sectionName='fiber', offset=0.0,
                               offsetType=MIDDLE_SURFACE, offsetField='',
                               thicknessAssignment=FROM_SECTION)

        # cohesive material
        model.CohesiveSection(name='cohesive',
                              material='material_cohesive', response=TRACTION_SEPARATION,
                              outOfPlaneThickness=None)
        part.SectionAssignment(region=part.sets['cohesive_face'], sectionName='cohesive', offset=0.0,
                               offsetType=MIDDLE_SURFACE, offsetField='',
                               thicknessAssignment=FROM_SECTION)

        sim_assembly = model.rootAssembly
        sim_instance = sim_assembly.Instance(
            dependent=ON, name=instance_name, part=part)

        # pbcs-----------------------------------------------------------------
        # right reference point
        right_reference_point_id = sim_assembly.ReferencePoint(
            point=(self.length * 1.5, self.center[1], 0.0)).id
        # top reference point
        top_reference_point_id = sim_assembly.ReferencePoint(
            point=(self.center[0], self.width * 1.5, 0.0)).id
        # reference points

        reference_points = sim_assembly.referencePoints
        # create sets for reference points
        sim_assembly.Set(
            name="Ref-R",
            referencePoints=((reference_points[right_reference_point_id],)),
        )
        sim_assembly.Set(
            name="Ref-T",
            referencePoints=((reference_points[top_reference_point_id],)),
        )

        # pbc for vertices ------------------------------------------------------------
        # regenerate assembly
        session.viewports["Viewport: 1"].setValues(
            displayedObject=sim_assembly)
        sim_assembly.regenerate()

        #
        set_name_for_vertices = ["NodeLB", "NodeLT", "NodeRB", "NodeRT"]
        for ii in range(len(vertices_name)):
            vertex_node = part.sets[vertices_name[ii]].nodes
            # create the assembly node set for vertex
            sim_assembly.SetFromNodeLabels(
                name=set_name_for_vertices[ii],
                nodeLabels=((instance_name, (vertex_node[0].label,)),),
                unsorted=True,
            )

        # apply the pbc to vertexes
        model.Equation(
            name="LB_RT_1",
            terms=(
                (1, "NodeRT", 1),
                (-1, "NodeLB", 1),
                (-1 * self.length, "Ref-R", 1),
                (-1 * self.width, "Ref-T", 1),
            ),
        )
        model.Equation(
            name="LB_RT_2",
            terms=(
                (1, "NodeRT", 2),
                (-1, "NodeLB", 2),
                (-1 * self.length, "Ref-R", 2),
                (-1 * self.width, "Ref-T", 2),
            ),
        )
        model.Equation(
            name="LT_RB_1",
            terms=(
                (1, "NodeRB", 1),
                (-1, "NodeLT", 1),
                (-1 * self.length, "Ref-R", 1),
                (1 * self.width, "Ref-T", 1),
            ),
        )
        model.Equation(
            name="LT_RB_2",
            terms=(
                (1, "NodeRB", 2),
                (-1, "NodeLT", 2),
                (-1 * self.length, "Ref-R", 2),
                (1 * self.width, "Ref-T", 2),
            ),
        )

        # pbcs for edges
        # part 1: equations for edges 2 (edgesFRONT_RIGHT) and 4 (edgesBACK_LEFT)
        edgesRIGHT_nodes = part.sets["edgesRIGHT"].nodes
        edgesRIGHT_nodes_sorted = sorted(edgesRIGHT_nodes, key=get_node_y)
        edgesLEFT_nodes = part.sets["edgesLEFT"].nodes
        edgesLEFT_nodes_sorted = sorted(edgesLEFT_nodes, key=get_node_y)
        if len(edgesRIGHT_nodes_sorted) == len(edgesLEFT_nodes_sorted):
            for ii in range(1, len(edgesRIGHT_nodes_sorted) - 1):
                sim_assembly.SetFromNodeLabels(
                    name="LEFT_" + str(ii),
                    nodeLabels=(
                        (
                            instance_name,
                            tuple([edgesLEFT_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                sim_assembly.SetFromNodeLabels(
                    name="RIGHT_" + str(ii),
                    nodeLabels=(
                        (
                            instance_name,
                            tuple([edgesRIGHT_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                for jj in range(1, 3):
                    model.Equation(
                        name="LEFT_RIGHT_" + str(ii) + "_" + str(jj),
                        terms=(
                            (1, "RIGHT_" + str(ii), jj),
                            (-1, "LEFT_" + str(ii), jj),
                            (-1 * self.length, "Ref-R", jj),
                        ),
                    )
        else:
            print "the number of nodes between the two sides are not the same"

        # part II:
        edgesTOP_nodes = part.sets["edgesTOP"].nodes
        edgesTOP_nodes_sorted = sorted(edgesTOP_nodes, key=get_node_x)
        edgesBOT_nodes = part.sets["edgesBOT"].nodes
        edgesBOT_nodes_sorted = sorted(edgesBOT_nodes, key=get_node_x)
        if len(edgesTOP_nodes_sorted) == len(edgesBOT_nodes_sorted):
            for ii in range(1, len(edgesTOP_nodes_sorted) - 1):
                sim_assembly.SetFromNodeLabels(
                    name="TOP_" + str(ii),
                    nodeLabels=(
                        (
                            instance_name,
                            tuple([edgesTOP_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                sim_assembly.SetFromNodeLabels(
                    name="BOT_" + str(ii),
                    nodeLabels=(
                        (
                            instance_name,
                            tuple([edgesBOT_nodes_sorted[ii].label]),
                        ),
                    ),
                    unsorted=True,
                )
                for jj in range(1, 3):
                    model.Equation(
                        name="TOP_BOT_" + str(ii) + "_" + str(jj),
                        terms=(
                            (1, "TOP_" + str(ii), jj),
                            (-1, "BOT_" + str(ii), jj),
                            (-1 * self.width, "Ref-T", jj),
                        ),
                    )
        else:
            print "the number of nodes between the two sides are not the same"

        # steps (static-step, implicit solver)
        model.StaticStep(name="Step-1", previous="Initial")
        step = model.StaticStep(
            initialInc=0.01,
            maxInc=1.0,
            maxNumInc=100000,
            minInc=1e-05,
            name="Step-1",
            previous="Initial",
            timePeriod=self.time_period,
            stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
            continueDampingFactors=True,
            adaptiveDampingRatio=0.0001
        )
        # # define filed output
        model.fieldOutputRequests["F-Output-1"].setValues(
            variables=(
                "S",
                "E",
                "LE",
                "ENER",
                "ELEN",
                "ELEDEN",
                "EVOL",
                "IVOL",
                'SDEG',
                'STATUS',
                'SDV'
            ),
            timeInterval=self.time_interval,
        )
        model.FieldOutputRequest(
            name="F-Output-2",
            createStepName="Step-1",
            variables=("U", "RF"),
            timeInterval=self.time_interval,
        )
        model.historyOutputRequests["H-Output-1"].setValues(
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

        # loadings --------------------------------------------------------------------
        model.DisplacementBC(
            name="E_11",
            createStepName="Step-1",
            region=sim_assembly.sets["Ref-R"],
            u1=self.strain[0],
            u2=UNSET,
            ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName="",
            localCsys=None,
        )
        # fix the rigid movement
        model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="rigid_x_1",
            region=sim_assembly.instances[instance_name].sets["VertexLB"],
            u1=0.0,
            u2=UNSET,
            ur3=UNSET,
        )

        model.DisplacementBC(
            amplitude=UNSET,
            createStepName="Step-1",
            distributionType=UNIFORM,
            fieldName="",
            fixed=OFF,
            localCsys=None,
            name="rigid_x_2",
            region=sim_assembly.instances[instance_name].sets["VertexLT"],
            u1=0.0,
            u2=0.0,
            ur3=UNSET,
        )

    def create_job(self):
        # run the job
        self.job = mdb.Job(name=self.job_name, model=self.model_name, description='', type=ANALYSIS,
                           atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
                           memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                           explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                           modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=self.subroutine_path,
                           scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=self.num_cpu,
                           numDomains=self.num_cpu, numGPUs=0)
        self.job.writeInput(consistencyChecking=OFF)

    def data_check(self):
        """check if there is error in the model
        """
        self.job.submit(consistencyChecking=OFF, datacheckJob=True)
        self.job.waitForCompletion()

    def submit_job(self):
        # submit the job
        self.job.submit(consistencyChecking=OFF)
        self.job.waitForCompletion()


class PostProcess:

    def __init__(self, dict):
        # job name
        self.job_name = str(dict["job_name"])
        # record time step (fir saving memory)
        self.record_time_step = dict["record_time_step"]
        # post process
        self.post_process()

    def post_process(self):
        # Define name of this .odb file
        odbfile = self.job_name + ".odb"
        # Open the output database
        rve_odb = openOdb(path=odbfile)
        # get element sets
        entire_element_set = rve_odb.rootAssembly.elementSets[" ALL ELEMENTS"]
        fiber_element_set = rve_odb.rootAssembly.instances["FINAL_STUFF"].elementSets['FIBERFACE']
        matrix_element_set = rve_odb.rootAssembly.instances["FINAL_STUFF"].elementSets["MATRIXFACE"]
        cohesive_element_set = rve_odb.rootAssembly.instances[
            "FINAL_STUFF"].elementSets["COHESIVE_FACE"]

        ref1_node_set = rve_odb.rootAssembly.nodeSets["REF-R"]
        ref2_node_set = rve_odb.rootAssembly.nodeSets["REF-T"]

        # get total steps, usually only one step
        my_steps = rve_odb.steps

        # get total frames
        total_frames = 0
        for ii in range(len(my_steps)):
            total_frames = total_frames + \
                len(my_steps[my_steps.keys()[ii]].frames)

        # get the variables that do not change with time steps
        rve_frame = rve_odb.steps[my_steps.keys()[0]].frames[0]

        # Extract volume at integration point in ENTIRE RVE:
        ivol_field = rve_frame.fieldOutputs["IVOL"]

        # get element volume for fiber
        self.ivol_fibers = self.get_ivol(ivol_field, fiber_element_set)
        # get the filed output for matrix
        self.ivol_matrix = self.get_ivol(ivol_field, matrix_element_set)
        # get the filed output for cohesive region
        self.ivol_cohesive = self.get_ivol(ivol_field, cohesive_element_set)

        # define required variables
        self.deformation_gradient = numpy.zeros((total_frames, 2, 2))
        self.strain = numpy.zeros((total_frames, 2, 2))
        self.stress = numpy.zeros((total_frames, 2, 2))
        self.total_time = numpy.zeros(len(my_steps))

        # define other variables
        self.U_ref1 = self.define_arrays(
            "U", rve_frame, ref1_node_set, total_frames)
        self.U_ref2 = self.define_arrays(
            "U", rve_frame, ref2_node_set, total_frames)
        self.RF_ref1 = self.define_arrays(
            "RF", rve_frame, ref1_node_set, total_frames)
        self.RF_ref2 = self.define_arrays(
            "RF", rve_frame, ref2_node_set, total_frames)
        self.plastic_strain = self.define_arrays(
            "SDV17", rve_frame, matrix_element_set, total_frames)

        # loop over all steps
        for ii in range(len(my_steps)):
            step = my_steps[my_steps.keys()[ii]]
            self.total_time[ii] = step.timePeriod
            # loop over all frames
            step_frames = step.frames  # improve speed
            for jj in range(len(step_frames)):
                frame = step_frames[jj]
                # get displacement for reference points
                # loop over elements
                u_field_ref1 = frame.fieldOutputs["U"].getSubset(
                    region=ref1_node_set, position=NODAL)
                u_field_ref2 = frame.fieldOutputs["U"].getSubset(
                    region=ref2_node_set, position=NODAL)
                if isinstance(u_field_ref1.values[0].data, float):
                    # variable is a scalar

                    self.U_ref1[ii * len(step_frames) +
                                jj] = u_field_ref1.values[0].data
                    self.U_ref2[ii * len(step_frames) +
                                jj] = u_field_ref2.values[0].data
                else:

                    # variable is an array
                    self.U_ref1[ii * len(step_frames) +
                                jj][:] = u_field_ref1.values[0].data[:]
                    self.U_ref2[ii * len(step_frames) +
                                jj][:] = u_field_ref2.values[0].data[:]

                # get reaction force for reference points
                rf_field_ref1 = frame.fieldOutputs["RF"].getSubset(
                    region=ref1_node_set, position=NODAL)
                rf_field_ref2 = frame.fieldOutputs["RF"].getSubset(
                    region=ref2_node_set, position=NODAL)
                if isinstance(rf_field_ref1.values[0].data, float):
                    # variable is a scalar
                    self.RF_ref1[ii * len(step_frames) +
                                 jj] = rf_field_ref1.values[0].data
                    self.RF_ref2[ii * len(step_frames) +
                                 jj] = rf_field_ref2.values[0].data
                else:

                    # variable is an array
                    self.RF_ref1[ii * len(step_frames) +
                                 jj][:] = rf_field_ref1.values[0].data[:]
                    self.RF_ref2[ii * len(step_frames) +
                                 jj][:] = rf_field_ref2.values[0].data[:]
                # get plastic strain for matrix

                # save the results every 10 frames to save memory
                if jj % self.record_time_step == 0:

                    plastic_strain_field = frame.fieldOutputs["SDV17"].getSubset(
                        region=matrix_element_set, position=INTEGRATION_POINT)
                    for kk in range(0, len(plastic_strain_field.values)):
                        self.plastic_strain[ii * len(step_frames) +
                                            jj / self.record_time_step][kk] = plastic_strain_field.values[kk].data
                # get deformation gradient
                for i in range(0, 2):
                    # get deformation gradient
                    self.deformation_gradient[ii * len(step_frames) + jj][0][i] \
                        = self.U_ref1[ii * len(step_frames) +
                                      jj][i] + numpy.identity(2)[0][i]
                    self.deformation_gradient[ii * len(step_frames) + jj][1][i] \
                        = self.U_ref2[ii * len(step_frames) +
                                      jj][i] + numpy.identity(2)[1][i]
                    # get strain
                    self.strain[ii * len(step_frames) + jj][0][i] \
                        = self.U_ref1[ii * len(step_frames) + jj][i]
                    self.strain[ii * len(step_frames) + jj][1][i] \
                        = self.U_ref2[ii * len(step_frames) + jj][i]

                    # get stress
                    self.stress[ii * len(step_frames) + jj][0][i] \
                        = self.RF_ref1[ii * len(step_frames) + jj][i] / \
                        (self.ivol_matrix.sum() + self.ivol_fibers.sum())
                    self.stress[ii * len(step_frames) + jj][1][i] = \
                        self.RF_ref2[ii * len(step_frames) + jj][i] / \
                        (self.ivol_matrix.sum() + self.ivol_fibers.sum())

        self.save_results()

    def get_ivol(self, field, element_set):

        # get the subfile
        ivolSubField = field.getSubset(
            region=element_set, position=INTEGRATION_POINT
        )
        # preallocate array for fibers
        ivol = numpy.zeros((len(ivolSubField.values)))
        # get the ivol
        for i in range(0, len(ivolSubField.values)):
            # Volume for i-th integration point
            ivol[i] = ivolSubField.values[i].data

        return ivol

    def define_arrays(self, field_name, frame, element_set, total_frames):
        # get the filed output
        field = frame.fieldOutputs[field_name]
        # preallocate array element set
        if field_name == "U" or field_name == "RF":

            # get filed output for selected element set
            sub_field = field.getSubset(region=element_set,
                                        position=NODAL)
            if isinstance(sub_field.values[0].data, float):
                # variable is a scalar
                array_temp = numpy.zeros((total_frames))
            else:
                # variable is an array
                array_temp = numpy.zeros((total_frames,
                                          len(sub_field.values[0].data)))
        elif field_name == "SDV17":
            # only save the results every 10 frames to save memory
            recorded_frames = int(
                ceil(total_frames / float(self.record_time_step)))
            sub_field = field.getSubset(region=element_set,
                                        position=INTEGRATION_POINT)
            array_temp = numpy.zeros((recorded_frames, len(sub_field.values)))

        return array_temp

    def save_results(self):

        # Save all variables to a single structured variable with all the data
        results = {
            "total_time": self.total_time,
            "stress": self.stress,
            "deformation_gradient": self.deformation_gradient,
            "strain": self.strain,
            "fiber_volume": self.ivol_fibers,
            "matrix_volume": self.ivol_matrix,
            "cohesive_volume": self.ivol_cohesive,
            "plastic_strain": self.plastic_strain,
        }
        # Save the results to a pickle file
        with open("results.pkl", "wb") as fp:
            pickle.dump(results, fp)


# get node coordinates
def get_node_y(node):
    return node.coordinates[1]

# get node coordinates


def get_node_x(node):
    return node.coordinates[0]
