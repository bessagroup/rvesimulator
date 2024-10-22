
# -*- coding: mbcs -*-
import assembly
import mesh
import numpy as np
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

class HyperelasticRVE(object):
    """ Hyperelastic RVE with hyperelastic matrix and hyperelastic inclusions
    """

    def __init__(self, sim_info={"inclusion_location_information": None,
                                 "radius_mu": None,
                                 "radius_std": None,
                                 "length_start": None,
                                 "length_end": None,
                                 "width_start": None,
                                 "width_end": None,
                                 "model_name": "rve_2phase",
                                 "part_name": "rve_2phase",
                                 "instance_name": "rve_2phase",
                                 "job_name": "hyperelastic_rve",
                                 "displacement_gradient": [[0.1, 0.0],[0.0, 0.0]],
                                 "params_matrix": None,
                                 "params_inclusion": None,
                                 "mesh_division": None,
                                 "simulation_time": 1.0,
                                 "num_pseudo_time_steps": None,
                                 "num_cpu": 1}):

        # ----------------------------parameters----------------------------- #

        # names of model, part, instance, job

        self.model_name = str(sim_info["model_name"])
        self.part_name = str(sim_info["part_name"])
        self.instance_name = str(sim_info["instance_name"])
        self.job_name = str(sim_info["job_name"])
        self.num_cpu = sim_info["num_cpu"]

        # information for creating geometry of the RVE

        self.inclusion_location_information = sim_info["inclusion_location_information"]
        self.length = (sim_info["length_end"] 
                       - sim_info["length_start"]
                       - 2 * sim_info["radius_mu"])
        self.width = (sim_info["width_end"]
                        - sim_info["width_start"]
                        - 2 * sim_info["radius_mu"])
        self.center = [ 
            (sim_info["length_end"] + sim_info["length_start"]) / 2.0,
            (sim_info["width_end"] + sim_info["width_start"]) / 2.0 ]
        self.radius_mu = sim_info["radius_mu"]
        self.radius_std = sim_info["radius_std"]

        # information for creating mesh 
        self.mesh_division = sim_info["mesh_division"]
        self.mesh_size = min(self.length, self.width) / self.mesh_division

        # information for simulation
        self.displacement_gradient = sim_info["displacement_gradient"]
        self.simulation_time = sim_info["simulation_time"]
        self.num_pseudo_time_steps = sim_info["num_pseudo_time_steps"]
        self.reporting_interval = self.simulation_time / self.num_pseudo_time_steps

        # information for material properties
        self.params_matrix = sim_info["params_matrix"]
        self.params_inclusion = sim_info["params_inclusion"]

        # submit # This is a strange way of implementing the code
        self.create_model()
        self.create_job()
        self.data_check()

    def create_model(self):
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

        # create the inclusions
        for ii in range(len(self.inclusion_location_information)):
            # actual inclusions
            sketch.CircleByCenterPerimeter(
                center=(self.inclusion_location_information[ii][0], self.inclusion_location_information[ii][1]),
                point1=(self.inclusion_location_information[ii][0] +
                        self.inclusion_location_information[ii][2], self.inclusion_location_information[ii][1]),
            )
        part.PartitionFaceBySketch(faces=faces[:], sketch=sketch)
        del sketch

        # create sets----------------------------------------------------------
        # get faces
        faces = part.faces[:]

        # for all faces

        part.Set(faces=faces, name="all_faces")
        inclusionface = part.faces.getByBoundingBox(
                    xMin=self.inclusion_location_information[0][0] - self.inclusion_location_information[0][2] -
                    0.001 * self.inclusion_location_information[0][2],
                    xMax=self.inclusion_location_information[0][0] + self.inclusion_location_information[0][2] +
                    0.001 * self.inclusion_location_information[0][2],
                    yMin=self.inclusion_location_information[0][1] - self.inclusion_location_information[0][2] -
                    0.001 * self.inclusion_location_information[0][2],
                    yMax=self.inclusion_location_information[0][1] + self.inclusion_location_information[0][2] +
                    0.001 * self.inclusion_location_information[0][2],
                    zMin=0.0,
                    zMax=1.0,
                )
        part.Set(faces=inclusionface, name="inclusionface")

        for ii in range(1, len(self.inclusion_location_information)):
            inclusionface_1 = part.faces.getByBoundingBox(
                xMin=self.inclusion_location_information[ii][0] - self.inclusion_location_information[ii][2] -
                0.001 * self.inclusion_location_information[ii][2],
                xMax=self.inclusion_location_information[ii][0] + self.inclusion_location_information[ii][2] +
                0.001 * self.inclusion_location_information[ii][2],
                yMin=self.inclusion_location_information[ii][1] - self.inclusion_location_information[ii][2] -
                0.001 * self.inclusion_location_information[ii][2],
                yMax=self.inclusion_location_information[ii][1] + self.inclusion_location_information[ii][2] +
                0.001 * self.inclusion_location_information[ii][2],
                zMin=0.0,
                zMax=1.0,
            )
            part.Set(faces=inclusionface_1, name="inclusionface_1")
            part.SetByBoolean(
                name="inclusionface",
                sets=(part.sets["inclusionface_1"],  part.sets["inclusionface"]),
                operation=UNION,
            )
        # delete fiber face 1
        del part.sets["inclusionface_1"]
        # create set for matrix
        part.SetByBoolean(
            name="matrixface",
            sets=(part.sets["all_faces"], part.sets["inclusionface"]),
            operation=DIFFERENCE,
        )

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
        part.seedPart(deviationFactor=0.05,  minSizeFactor=0.05, size=self.mesh_size)
        # element type (plane strain element)
        elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD,
                                secondOrderAccuracy=OFF, hourglassControl=ENHANCED,
                                distortionControl=DEFAULT)
        elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)
        part.setElementType(
            regions=part.sets["inclusionface"], elemTypes=(elemType1, elemType2))
        part.setElementType(
            regions=part.sets["matrixface"], elemTypes=(elemType1, elemType2))

        part.generateMesh()

        # material properties ---------------------------------------------------------
        # material property for particles
        material_inclusion = model.Material(name="material_inclusion")
        material_inclusion.Hyperelastic(materialType=ISOTROPIC, testData=OFF, type=NEO_HOOKE, 
                                        volumetricResponse=VOLUMETRIC_DATA, 
                                        table=((self.params_inclusion[0], self.params_inclusion[1]), ))

        # material property for the matrix part
        material_matrix = model.Material(name="material_matrix")
        material_matrix.Hyperelastic(materialType=ISOTROPIC, testData=OFF, type=ARRUDA_BOYCE, 
                                        volumetricResponse=VOLUMETRIC_DATA, 
                                        table=((self.params_matrix[0], self.params_matrix[1], self.params_matrix[2]), ))

        # create section and assign material property to corresponding section
        # matrix material
        model.HomogeneousSolidSection(
            name='matrix', material='material_matrix', thickness=None)
        part.SectionAssignment(region=part.sets['matrixface'], sectionName='matrix', offset=0.0,
                            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


        # inclusion material
        model.HomogeneousSolidSection(
            name='inclusion', material='material_inclusion', thickness=None)
        part.SectionAssignment(region=part.sets['inclusionface'], sectionName='inclusion', offset=0.0,
                            offsetType=MIDDLE_SURFACE, offsetField='',
                            thicknessAssignment=FROM_SECTION)

        sim_assembly = model.rootAssembly
        sim_instance = sim_assembly.Instance(
            dependent=ON, name=instance_name, part=part)
        
        # pbcs-----------------------------------------------------------------
        # right reference point
        right_reference_point_id = sim_assembly.ReferencePoint(point=(self.length * 1.5, self.center[1], 0.0)).id
        # top reference point
        top_reference_point_id = sim_assembly.ReferencePoint(point=(self.center[0], self.width * 1.5, 0.0)).id
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
        session.viewports["Viewport: 1"].setValues( displayedObject=sim_assembly)
        sim_assembly.regenerate()

        #
        set_name_for_vertices =["NodeLB", "NodeLT", "NodeRB", "NodeRT"]
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
           raise ValueError(
                "the number of nodes between the left and right are not the same"
            )

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
            raise ValueError(
                "the number of nodes between top and the bottom is not the same"
            )

        
        # steps (static-step, implicit solver) 
        model.StaticStep(name="Step-1", previous="Initial")
        step = model.StaticStep(
                    initialInc=0.01,
                    maxInc=1.0,
                    maxNumInc=100000,
                    minInc=1e-20,
                    name="Step-1",
                    previous="Initial",
                    timePeriod=self.simulation_time,
                    nlgeom=ON
                )
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
                    ),
                    timeInterval=self.reporting_interval,
                )
        model.FieldOutputRequest(
                    name="F-Output-2",
                    createStepName="Step-1",
                    variables=("U", "RF"),
                    timeInterval=self.reporting_interval,
                )


        # loadings --------------------------------------------------------------------
        model.DisplacementBC(
                    name="RightBC",
                    createStepName="Step-1",
                    region=sim_assembly.sets["Ref-R"],
                    u1=self.displacement_gradient[0][0],
                    u2=self.displacement_gradient[0][1],
                    ur3=UNSET,
                    amplitude=UNSET,
                    fixed=OFF,
                    distributionType=UNIFORM,
                    fieldName="",
                    localCsys=None,
                )
        
        model.DisplacementBC(
                    name="TopBC",
                    createStepName="Step-1",
                    region=sim_assembly.sets["Ref-T"],
                    u1=self.displacement_gradient[1][0],
                    u2=self.displacement_gradient[1][1],
                    ur3=UNSET,
                    amplitude=UNSET,
                    fixed=OFF,
                    distributionType=UNIFORM,
                    fieldName="",
                    localCsys=None,
                )
        
    def create_job(self):
        # run the job
        self.job = mdb.Job(name=self.job_name, model=self.model_name, description='', type=ANALYSIS,
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
            scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=self.num_cpu,
            numDomains=self.num_cpu, numGPUs=0)
        
    def data_check(self):
        """check if there is error in the model
        """
        self.job.submit(consistencyChecking=OFF, datacheckJob=True)  
        self.job.waitForCompletion()




class PostProcess:

    def __init__(self, dict):
        # job name
        job_name = dict["job_name"]
        self.job_name = str(job_name)
        self.post_process()

    def post_process(self):
        # Define name of this .odb file
        odbfile = self.job_name + ".odb"
        # Open the output database
        rve_odb = openOdb(path=odbfile)
        # get element sets
        entire_element_set = rve_odb.rootAssembly.elementSets[" ALL ELEMENTS"]
        inclusion_element_set = rve_odb.rootAssembly.instances["RVE_2PHASE"].elementSets['INCLUSIONFACE']
        matrix_element_set = rve_odb.rootAssembly.instances["RVE_2PHASE"].elementSets["MATRIXFACE"]

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

        # get element volume for inclusion
        self.ivol_inclusions = self.get_ivol(ivol_field, inclusion_element_set)
        # get the filed output for matrix
        self.ivol_matrix = self.get_ivol(ivol_field, matrix_element_set)

        # define required variables
        self.deformation_gradient = np.zeros((total_frames, 2, 2))
        self.displacement_gradient = np.zeros((total_frames, 2, 2))
        self.pk1_stress = np.zeros((total_frames, 2, 2))
        self.total_time = np.zeros(len(my_steps))

        # define other variables
        self.U_ref1 = self.define_arrays(
            "U", rve_frame, ref1_node_set, total_frames)
        self.U_ref2 = self.define_arrays(
            "U", rve_frame, ref2_node_set, total_frames)
        self.RF_ref1 = self.define_arrays(
            "RF", rve_frame, ref1_node_set, total_frames)
        self.RF_ref2 = self.define_arrays(
            "RF", rve_frame, ref2_node_set, total_frames)

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
                
                # get deformation gradient, displacement_gradient and pk1_stress
                for i in range(0, 2):
                    # get deformation gradient
                    self.deformation_gradient[ii * len(step_frames) + jj][0][i] \
                        = self.U_ref1[ii * len(step_frames) + jj][i] \
                            + np.identity(2)[0][i]
                    self.deformation_gradient[ii * len(step_frames) + jj][1][i] \
                        = self.U_ref2[ii * len(step_frames) + jj][i] \
                            + np.identity(2)[1][i]
                    # get displacement_gradient
                    self.displacement_gradient[ii * len(step_frames) + jj][0][i] \
                        = self.U_ref1[ii * len(step_frames) + jj][i]
                    self.displacement_gradient[ii * len(step_frames) + jj][1][i] \
                        = self.U_ref2[ii * len(step_frames) + jj][i]

                    # get pk1_stress
                    self.pk1_stress[ii * len(step_frames) + jj][0][i] \
                        = self.RF_ref1[ii * len(step_frames) + jj][i] / \
                        (self.ivol_matrix.sum() + self.ivol_inclusions.sum())
                    self.pk1_stress[ii * len(step_frames) + jj][1][i] = \
                        self.RF_ref2[ii * len(step_frames) + jj][i] / \
                        (self.ivol_matrix.sum() + self.ivol_inclusions.sum())

        self.save_results()

    def get_ivol(self, field, element_set):

        # get the subfile
        ivolSubField = field.getSubset(
            region=element_set, position=INTEGRATION_POINT
        )
        # preallocate array for inclusions
        ivol = np.zeros((len(ivolSubField.values)))
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
                array_temp = np.zeros((total_frames))
            else:
                # variable is an array
                array_temp = np.zeros((total_frames,
                                          len(sub_field.values[0].data)))

        return array_temp

    def save_results(self):

        # Save all variables to a single structured variable with all the data
        results = {
            "total_time": self.total_time,
            "pk1_stress": self.pk1_stress,
            "deformation_gradient": self.deformation_gradient,
            "displacement_gradient": self.displacement_gradient,
            "inclusion_volume": self.ivol_inclusions,
            "matrix_volume": self.ivol_matrix,
        }
        # Save the results to a pickle file
        with open("results.p", "wb") as fp:
            pickle.dump(results, fp)


# get node coordinates
def get_node_y(node):
    return node.coordinates[1]

# get node coordinates
def get_node_x(node):
    return node.coordinates[0]
