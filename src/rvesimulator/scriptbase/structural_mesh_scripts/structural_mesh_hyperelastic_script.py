# -*- coding: mbcs -*-
import assembly
import mesh
import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *

try:
    import cPickle as pickle  # Improve speed
except:
    import pickle


class HyperelasticRVE(object):
    """ Hyperelastic RVE with hyperelastic matrix and hyperelastic inclusions
    with structural mesh
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
                                 "displacement_gradient": [[0.1, 0.0], [0.0, 0.0]],
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
            (sim_info["width_end"] + sim_info["width_start"]) / 2.0]
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
    # ==========================  Meshing ========================== #
        part.seedPart(size=self.mesh_size,
                      deviationFactor=0.1, minSizeFactor=0.1)
        # Apply mesh controls for structured meshing using plane strain elements
        elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
        # Assign element type to the part
        part.setElementType(regions=(part.faces,),
                            elemTypes=(elemType1, elemType2))
        # Generate the mesh
        part.generateMesh()

        # create sets----------------------------------------------------------
        # get faces
        faces = part.faces[:]
        # for all faces
        part.Set(faces=faces, name="all_faces")

        # load the microstructure descriptor from a file
        microstructure_descriptor = np.load('microstructure.rgmsh.npy')
        microstructure_descriptor = microstructure_descriptor.flatten()
        # get the index with the value 0
        index = np.where(microstructure_descriptor == 1)
        index = index[0] + 1
        # create the set for the elements with the value 0
        part.Set(name='matrix',
                 elements=part.elements.sequenceFromLabels(index.tolist()))

        # get the index with the value 1
        index = np.where(microstructure_descriptor == 2)
        index = index[0] + 1
        # create the set for the elements with the value 1
        part.Set(name='fiber',
                 elements=part.elements.sequenceFromLabels(index.tolist()))

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

        # material properties ---------------------------------------------------------
        # material property for particles
        material_inclusion = model.Material(name="fiber")
        material_inclusion.Hyperelastic(materialType=ISOTROPIC,
                                        testData=OFF,
                                        type=NEO_HOOKE,
                                        volumetricResponse=VOLUMETRIC_DATA,
                                        table=((self.params_inclusion[0],
                                                self.params_inclusion[1]), ))
        # create section and assign material property to corresponding section
        model.HomogeneousSolidSection(
            name='fiber', material='fiber', thickness=None)
        part.SectionAssignment(region=part.sets['fiber'],
                               sectionName='fiber', offset=0.0,
                               offsetType=MIDDLE_SURFACE, offsetField='',
                               thicknessAssignment=FROM_SECTION)

        # material property for the matrix part
        material_matrix = model.Material(name="matrix")
        material_matrix.Hyperelastic(materialType=ISOTROPIC,
                                     testData=OFF,
                                     type=ARRUDA_BOYCE,
                                     volumetricResponse=VOLUMETRIC_DATA,
                                     table=((self.params_matrix[0],
                                             self.params_matrix[1],
                                             self.params_matrix[2]), ))
        # create section and assign material property to corresponding section
        model.HomogeneousSolidSection(
            name='matrix', material='matrix', thickness=None)
        part.SectionAssignment(region=part.sets['matrix'],
                               sectionName='matrix', offset=0.0,
                               offsetType=MIDDLE_SURFACE,
                               offsetField='',
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

        # pbc for vertices ----------------------------------------------------
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
        # pop the first and last node
        edgesRIGHT_nodes_sorted.pop(0)
        edgesRIGHT_nodes_sorted.pop(-1)

        # sort the nodes of the right edge
        edgesLEFT_nodes = part.sets["edgesLEFT"].nodes
        edgesLEFT_nodes_sorted = sorted(edgesLEFT_nodes, key=get_node_y)
        # pop the first and last node
        edgesLEFT_nodes_sorted.pop(0)
        edgesLEFT_nodes_sorted.pop(-1)

        if len(edgesRIGHT_nodes_sorted) != len(edgesLEFT_nodes_sorted):

            print("The number of nodes between the two sides are not the same")
            print(
                f"Right side node is: {len(edgesRIGHT_nodes_sorted)} ")
            print(
                f"Left side node is: {len(edgesLEFT_nodes_sorted)} ")
            if len(edgesRIGHT_nodes_sorted) > len(edgesLEFT_nodes_sorted):
                print(
                    f"remove {np.abs(len(edgesRIGHT_nodes_sorted) - len(edgesLEFT_nodes_sorted))} nodes from right ")
                num_remove = len(edgesRIGHT_nodes_sorted) - \
                    len(edgesLEFT_nodes_sorted)
                index_remove = np.random.choice(
                    len(edgesRIGHT_nodes_sorted) - num_remove, num_remove, replace=False)
                for jj in index_remove:
                    edgesRIGHT_nodes_sorted.pop(jj)
            else:
                print(
                    f"remove {np.abs(len(edgesRIGHT_nodes_sorted) - len(edgesLEFT_nodes_sorted))} nodes from left ")
                num_remove = len(edgesLEFT_nodes_sorted) - \
                    len(edgesRIGHT_nodes_sorted)
                index_remove = np.random.choice(
                    len(edgesLEFT_nodes_sorted) - num_remove, num_remove, replace=False)
                for jj in index_remove:
                    edgesLEFT_nodes_sorted.pop(jj)
        else:
            print("The number of nodes between the two sides are the same")

        # pbcs for the left and right edges
        for ii in range(len(edgesRIGHT_nodes_sorted)):
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
        # part II: pbcs for the top and bottom edges
        # get the top nodes
        edgesTOP_nodes = part.sets["edgesTOP"].nodes
        edgesTOP_nodes_sorted = sorted(edgesTOP_nodes, key=get_node_x)
        # pop the first and last node
        edgesTOP_nodes_sorted.pop(0)
        edgesTOP_nodes_sorted.pop(-1)

        # get the bottom nodes
        edgesBOT_nodes = part.sets["edgesBOT"].nodes
        edgesBOT_nodes_sorted = sorted(edgesBOT_nodes, key=get_node_x)
        # pop the first and last node
        edgesBOT_nodes_sorted.pop(0)
        edgesBOT_nodes_sorted.pop(-1)

        if len(edgesTOP_nodes_sorted) != len(edgesBOT_nodes_sorted):
            print("The number of nodes between the top and the bottom are not the same")
            print(
                f"Top side node is: {len(edgesTOP_nodes_sorted)} ")
            print(
                f"Bottom side node is: {len(edgesBOT_nodes_sorted)} ")
            if len(edgesTOP_nodes_sorted) > len(edgesBOT_nodes_sorted):
                print(
                    f"remove {np.abs(len(edgesTOP_nodes_sorted) - len(edgesBOT_nodes_sorted))} nodes from top ")
                num_remove = len(edgesTOP_nodes_sorted) - \
                    len(edgesBOT_nodes_sorted)
                index_remove = np.random.choice(
                    len(edgesTOP_nodes_sorted) - num_remove, num_remove, replace=False)
                for jj in index_remove:
                    edgesTOP_nodes_sorted.pop(jj)
            else:
                print(
                    f"remove {np.abs(len(edgesTOP_nodes_sorted) - len(edgesBOT_nodes_sorted))} nodes from bottom ")
                num_remove = len(edgesBOT_nodes_sorted) - \
                    len(edgesTOP_nodes_sorted)
                index_remove = np.random.choice(
                    len(edgesBOT_nodes_sorted) - num_remove, num_remove, replace=False)
                for jj in index_remove:
                    edgesBOT_nodes_sorted.pop(jj)
        else:
            print("The number of nodes between the top and the bottom are the same")

        # pbcs for the top and bottom edges
        for ii in range(len(edgesTOP_nodes_sorted)):
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

        # steps (static-step, implicit solver)
        model.StaticStep(name="Step-1", previous="Initial")
        step = model.StaticStep(
            initialInc=0.01,
            maxInc=1.0,
            maxNumInc=100000,
            minInc=1e-10,
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
        inclusion_element_set = rve_odb.rootAssembly.instances[
            "RVE_2PHASE"].elementSets['FIBER']
        matrix_element_set = rve_odb.rootAssembly.instances["RVE_2PHASE"].elementSets["MATRIX"]

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
        with open("results.pkl", "wb") as fp:
            pickle.dump(results, fp)


# get node coordinates
def get_node_y(node):
    return node.coordinates[1]


def get_node_x(node):
    return node.coordinates[0]
