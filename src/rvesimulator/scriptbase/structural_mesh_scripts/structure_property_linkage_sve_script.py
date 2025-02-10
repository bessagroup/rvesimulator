# -*- coding: mbcs -*-
import numpy as np
import pickle
from abaqus import *
from abaqusConstants import *
from caeModules import *
# import packages for abaqus post-processing
from odbAccess import *


def simulation_script(sim_info):
    # Create a model
    # names of model, part, instance
    model_name = "rve"
    part_name = "Final_Stuff"
    instance_name = "Final_Stuff"
    job_name = str(sim_info["job_name"])
    num_cpu = sim_info["num_cpu"]

    length = (
        sim_info["len_end"] - sim_info["len_start"]
    ) - 2 * sim_info["radius_mu"]
    width = (
        sim_info["wid_end"] - sim_info["wid_start"]
    ) - 2 * sim_info["radius_mu"]
    center = [
        (sim_info["len_end"] + sim_info["len_start"]) / 2.0,
        (sim_info["wid_end"] + sim_info["wid_start"]) / 2.0,
    ]

    # material properties (elastic fibers)
    youngs_modulus_fiber = sim_info["youngs_modulus_fiber"]
    poisson_ratio_fiber = sim_info["poisson_ratio_fiber"]
    # material properties (plastic matrix)
    youngs_modulus_matrix = sim_info["youngs_modulus_matrix"]
    poisson_ratio_matrix = sim_info["poisson_ratio_matrix"]
    hardening_table_matrix = np.zeros(
        (len(sim_info["hardening_table_matrix"][0]), 2))
    for ii in range(2):
        hardening_table_matrix[:, ii] = sim_info["hardening_table_matrix"][ii]

    # information of RVE modeling
    strain = sim_info["strain"]
    mesh_size = (
        min(length, width) / sim_info["mesh_partition"]
    )
    time_period = sim_info["simulation_time"]
    time_interval = (
        sim_info["simulation_time"] / sim_info["num_steps"]
    )
    delta = min(min(length, width) / 1000, mesh_size / 10)

    # create model---------------------------------------------------------
    mdb.models.changeKey(fromName='Model-1', toName=model_name)
    model = mdb.models[model_name]
    # delete existed model
    if 'Model-1' in mdb.models.keys():
        del mdb.models['Model-1']

    # create sketch--------------------------------------------------------
    # model the matrix part first
    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    sketch.rectangle(point1=(0.0, 0.0), point2=(length, width))
    part = model.Part(name=part_name,
                      dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sketch)
    del sketch

    # ==========================  Meshing ========================== #
    part.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    # Apply mesh controls for structured meshing
    part.setMeshControls(regions=part.faces,
                         elemShape=QUAD, technique=STRUCTURED)
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
    part.Set(name='matrix', elements=part.elements.sequenceFromLabels(index.tolist()))

    # get the index with the value 1
    index = np.where(microstructure_descriptor == 2)
    index = index[0] + 1
    # create the set for the elements with the value 1
    part.Set(name='fiber', elements=part.elements.sequenceFromLabels(index.tolist()))

    # create sets for edges
    s = part.edges
    edgesLEFT = s.getByBoundingBox(
        center[0] - length / 2 - delta,
        center[1] - width / 2 - delta,
        0,
        center[0] - length / 2 + delta,
        center[1] + width / 2 + delta,
        0,
    )
    part.Set(edges=edgesLEFT, name="edgesLEFT")
    edgesRIGHT = s.getByBoundingBox(
        center[0] + length / 2 - delta,
        center[1] - width / 2 - delta,
        0,
        center[0] + length / 2 + delta,
        center[1] + width / 2 + delta,
        0,
    )
    part.Set(edges=edgesRIGHT, name="edgesRIGHT")
    edgesTOP = s.getByBoundingBox(
        center[0] - length / 2 - delta,
        center[1] + width / 2 - delta,
        0,
        center[0] + width / 2 + delta,
        center[1] + width / 2 + delta,
        0,
    )
    part.Set(edges=edgesTOP, name="edgesTOP")
    edgesBOT = s.getByBoundingBox(
        center[0] - length / 2 - delta,
        center[1] - width / 2 - delta,
        0,
        center[0] + length / 2 + delta,
        center[1] - width / 2 + delta,
        0,
    )
    part.Set(edges=edgesBOT, name="edgesBOT")
    name_edges = ["edgesLEFT", "edgesRIGHT", "edgesTOP", "edgesBOT"]

    # create set for vertices
    v = part.vertices
    vertexLB = v.getByBoundingBox(
        center[0] - length / 2 - delta,
        center[1] - width / 2 - delta,
        0,
        center[0] - length / 2 + delta,
        center[1] - width / 2 + delta,
        0,
    )
    part.Set(vertices=vertexLB, name="VertexLB")
    vertexRB = v.getByBoundingBox(
        center[0] + length / 2 - delta,
        center[1] - width / 2 - delta,
        0,
        center[0] + length / 2 + delta,
        center[1] - width / 2 + delta,
        0,
    )
    part.Set(vertices=vertexRB, name="VertexRB")
    vertexRT = v.getByBoundingBox(
        center[0] + length / 2 - delta,
        center[1] + width / 2 - delta,
        0,
        center[0] + length / 2 + delta,
        center[1] + width / 2 + delta,
        0,
    )
    part.Set(vertices=vertexRT, name="VertexRT")
    vertexLT = v.getByBoundingBox(
        center[0] - length / 2 - delta,
        center[1] + width / 2 - delta,
        0,
        center[0] - length / 2 + delta,
        center[1] + width / 2 + delta,
        0,
    )
    part.Set(vertices=vertexLT, name="VertexLT")
    vertices_name = ["VertexLB", "VertexLT", "VertexRB", "VertexRT"]

    # material properties ---------------------------------------------------------
    # material property for fiber part
    material_fiber = model.Material(name="fiber")
    material_fiber.Elastic(
        table=((youngs_modulus_fiber, poisson_ratio_fiber),))
    model.HomogeneousSolidSection(
        name="fiber", material="fiber", thickness=None
    )
    part.SectionAssignment(
        region=part.sets["fiber"],
        sectionName="fiber",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )
    # material property for matrix part
    material_matrix = model.Material(name="matrix")
    material_matrix.Elastic(
        table=((youngs_modulus_matrix, poisson_ratio_matrix),))
    material_matrix.Plastic(table=(hardening_table_matrix))
    model.HomogeneousSolidSection(
        name="matrix", material="matrix", thickness=None
    )
    part.SectionAssignment(
        region=part.sets["matrix"],
        sectionName="matrix",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )

    sim_assembly = model.rootAssembly
    sim_instance = sim_assembly.Instance(
        dependent=ON, name=instance_name, part=part)

    # pbcs-----------------------------------------------------------------
    # right reference point
    right_reference_point_id = sim_assembly.ReferencePoint(
        point=(length * 1.5,  center[1], 0.0)).id
    # top reference point
    top_reference_point_id = sim_assembly.ReferencePoint(
        point=(center[0],  width * 1.5, 0.0)).id
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
            (-1 * length, "Ref-R", 1),
            (-1 * width, "Ref-T", 1),
        ),
    )
    model.Equation(
        name="LB_RT_2",
        terms=(
            (1, "NodeRT", 2),
            (-1, "NodeLB", 2),
            (-1 * length, "Ref-R", 2),
            (-1 * width, "Ref-T", 2),
        ),
    )
    model.Equation(
        name="LT_RB_1",
        terms=(
            (1, "NodeRB", 1),
            (-1, "NodeLT", 1),
            (-1 * length, "Ref-R", 1),
            (1 * width, "Ref-T", 1),
        ),
    )
    model.Equation(
        name="LT_RB_2",
        terms=(
            (1, "NodeRB", 2),
            (-1, "NodeLT", 2),
            (-1 * length, "Ref-R", 2),
            (1 * width, "Ref-T", 2),
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
                        (-1 * length, "Ref-R", jj),
                    ),
                )
    else:
        raise ValueError(
            "the number of nodes between the two sides are not the same")

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
                        (-1 * width, "Ref-T", jj),
                    ),
                )
    else:
        raise ValueError(
            "the number of nodes between the two sides are not the same")

    #  ==========================  Create step ========================== #
    model.StaticStep(name="Step-1", previous="Initial")
    step = model.StaticStep(
        initialInc=0.01,
        maxInc=1.0,
        maxNumInc=100000,
        minInc=1e-20,
        name="Step-1",
        previous="Initial",
        timePeriod=time_period,
        nlgeom=ON,
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
            "PEEQ",

        ),
        timeInterval=time_interval,
    )
    model.FieldOutputRequest(
        name="F-Output-2",
        createStepName="Step-1",
        variables=("U", "RF"),
        timeInterval=time_interval,
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
        timeInterval=time_interval,
    )

    # ==========================  Apply boundary conditions ========================== #

    if strain[1] == 0 and strain[2] == 0:
        # only E11
        E11(strain[0], model, sim_assembly)
        # rigid body motion for x direction
        rigid_constraints_for_x_direction(model, sim_assembly, instance_name)
    elif strain[0] == 0 and strain[2] == 0:
        # only E22
        E22(E22, model, sim_assembly)
        # rigid body motion for y direction
        rigid_constraints_for_y_direction(model, sim_assembly, instance_name)
    elif strain[0] == 0 and strain[1] == 0:
        # only E12
        E12(strain[2], model, sim_assembly)
    else:
        E11(strain[0], model, sim_assembly)
        E22(strain[1], model, sim_assembly)
        E12(strain[2], model, sim_assembly)

    # ==========================  Create job ========================== #
    # create job
    mdb.Job(
        name=job_name,
        model=model_name,
        description="",
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine="",
        scratch="",
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=num_cpu,
        numGPUs=0,
        numDomains=num_cpu)
    # submit the job
    mdb.jobs[job_name].writeInput(consistencyChecking=OFF)


def E11(E11, model, assembly):
    "displacement for x direction"
    model.DisplacementBC(
        name="E_11",
        createStepName="Step-1",
        region=assembly.sets["Ref-R"],
        u1=E11,
        u2=UNSET,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )


def E22(E22, model, assembly):
    "displacement for y direction"
    model.DisplacementBC(
        name="E_22",
        createStepName="Step-1",
        region=assembly.sets["Ref-T"],
        u1=UNSET,
        u2=E22,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )


def E12(E12, model, assembly):
    "displacement for xy direction"
    model.DisplacementBC(
        name="E_12",
        createStepName="Step-1",
        region=assembly.sets["Ref-R"],
        u1=UNSET,
        u2=E12,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )
    model.DisplacementBC(
        name="E_21",
        createStepName="Step-1",
        region=assembly.sets["Ref-T"],
        u1=E12,
        u2=UNSET,
        ur3=UNSET,
        amplitude=UNSET,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName="",
        localCsys=None,
    )


def rigid_constraints_for_x_direction(model, assembly, instance_name):
    """displacement rigid constraints for x direction"""
    model.DisplacementBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="rigid_x_1",
        region=assembly.instances[instance_name].sets["VertexLB"],
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
        region=assembly.instances[instance_name].sets["VertexLT"],
        u1=0.0,
        u2=0.0,
        ur3=UNSET,
    )


def rigid_constraints_for_y_direction(model, assembly, instance_name):
    model.DisplacementBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="rigid_y_1",
        region=assembly.instances[instance_name].sets["VertexLB"],
        u1=0.0,
        u2=0.0,
        ur3=UNSET,
    )
    model.DisplacementBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="rigid_y_2",
        region=assembly.instances[instance_name].sets["VertexRB"],
        u1=UNSET,
        u2=0.0,
        ur3=UNSET,
    )

# get node coordinates


def get_node_y(node):
    return node.coordinates[1]


def get_node_x(node):
    return node.coordinates[0]


def post_process(dict):
    """get the rve results"""

    # Define the name of the Phase_1 part
    odbfile = str(dict["job_name"]) + ".odb"  # Define name of this .odb file
    RVEodb = openOdb(path=odbfile)
    # Determine the number of steps in the output database.
    mySteps = RVEodb.steps
    numSteps = len(mySteps)

    entireRVE_elSet = RVEodb.rootAssembly.elementSets[" ALL ELEMENTS"]

    dummy1_nSet = RVEodb.rootAssembly.nodeSets["REF-R"]
    dummy2_nSet = RVEodb.rootAssembly.nodeSets["REF-T"]
    #
    # For each step, obtain the following:
    #     1) The step key.
    #     2) The number of frames in the step.
    #     3) The increment number of the last frame in the step.
    #
    totalNumFrames = 0
    for iStep in range(numSteps):
        stepKey = mySteps.keys()[iStep]
        step = mySteps[stepKey]
        numFrames = len(step.frames)
        totalNumFrames = totalNumFrames + numFrames
    # getting the history output
    ALLPD_data = (
        RVEodb.steps[mySteps.keys()[0]]
        .historyRegions["Assembly ASSEMBLY"]
        .historyOutputs["ALLPD"]
        .data
    )
    ALLPD = np.zeros((len(ALLPD_data), 2))
    for ii in range(len(ALLPD_data)):
        ALLPD[ii, 0] = ALLPD_data[ii][0]
        ALLPD[ii, 1] = ALLPD_data[ii][1]

    # Preallocate quantities for speed
    RVEframe = RVEodb.steps[mySteps.keys()[0]].frames[
        0
    ]  # Undeformed config.
    # Extract volume at integration point in ENTIRE RVE:
    ivolField = RVEframe.fieldOutputs["IVOL"]
    ivolSubField = ivolField.getSubset(
        region=entireRVE_elSet, position=INTEGRATION_POINT
    )
    ivol = np.zeros((len(ivolSubField.values)))
    # stress = numpy.zeros((len(ivolSubField.values)))
    tot_vol = 0.0
    for i in range(0, len(ivolSubField.values)):
        ivol[i] = ivolSubField.values[
            i
        ].data  # Volume for i-th integration point
        tot_vol = tot_vol + ivol[i]  # total volume

    # define arrays to store the displacement ========================== #
    U_Field = RVEframe.fieldOutputs["U"]
    U_dummy1_SubField = U_Field.getSubset(
        region=dummy1_nSet, position=NODAL
    )
    U_dummy2_SubField = U_Field.getSubset(
        region=dummy2_nSet, position=NODAL
    )
    #
    if isinstance(U_dummy1_SubField.values[0].data, float):
        # Then variable is a scalar
        U_dummy1 = np.zeros((totalNumFrames))
        U_dummy2 = np.zeros((totalNumFrames))
    else:
        # Variable is an array
        U_dummy1 = np.zeros(
            (totalNumFrames, len(U_dummy1_SubField.values[0].data))
        )
        U_dummy2 = np.zeros(
            (totalNumFrames, len(U_dummy2_SubField.values[0].data))
        )
    # define arrays to store the reaction forces ====================== #
    RF_Field = RVEframe.fieldOutputs["RF"]
    RF_dummy1_SubField = RF_Field.getSubset(
        region=dummy1_nSet, position=NODAL
    )
    RF_dummy2_SubField = RF_Field.getSubset(
        region=dummy2_nSet, position=NODAL
    )

    #
    if isinstance(RF_dummy1_SubField.values[0].data, float):
        # Then variable is a scalar
        RF_dummy1 = np.zeros((totalNumFrames))
        RF_dummy2 = np.zeros((totalNumFrames))
    else:
        # Variable is an array
        RF_dummy1 = np.zeros(
            (totalNumFrames, len(RF_dummy1_SubField.values[0].data))
        )
        RF_dummy2 = np.zeros(
            (totalNumFrames, len(RF_dummy2_SubField.values[0].data))
        )

    # define array to save Von Mises stress ============================ #
    S_Field = RVEframe.fieldOutputs["S"]
    S_SubField = S_Field.getSubset(
        region=entireRVE_elSet, position=INTEGRATION_POINT
    )
    #
    if isinstance(S_SubField.values[0].mises, float):
        # Then variable is a scalar
        VonMises = np.zeros((totalNumFrames))
    else:
        # Variable is an array
        VonMises = np.zeros(
            (totalNumFrames, len(S_SubField.values[0].mises))
        )

    # define array to save Plastic Energy ============================ #
    PEEQ_Field = RVEframe.fieldOutputs["PEEQ"]
    PEEQ_SubField = PEEQ_Field.getSubset(
        region=entireRVE_elSet, position=INTEGRATION_POINT
    )
    #
    if isinstance(PEEQ_SubField.values[0].data, float):
        # Then variable is a scalar
        PEEQ = np.zeros((totalNumFrames))
    else:
        # Variable is an array
        PEEQ = np.zeros(
            (totalNumFrames, len(PEEQ_SubField.values[0].data))
        )

    # Loop over Steps and Frames to compute average quantities in RVE
    eye = np.identity(2)
    deformation_gradient = np.zeros((totalNumFrames, 2, 2))
    strain = np.zeros((totalNumFrames, 2, 2))
    nominal_stress = np.zeros((totalNumFrames, 2, 2))
    total_time = np.zeros(numSteps)
    previousFrame = 0
    numFrames = 0
    for iStep in range(numSteps):
        previousFrame = previousFrame + numFrames
        stepKey = mySteps.keys()[iStep]
        step = mySteps[stepKey]
        total_time[iStep] = step.timePeriod
        numFrames = len(step.frames)
        #
        for iFrame_step in range(0, numFrames):
            iFrame = previousFrame + iFrame_step
            RVEframe = RVEodb.steps[stepKey].frames[iFrame]
            # Finished computing average for this variable!
            # Variable: U
            U_Field = RVEframe.fieldOutputs["U"]
            U_dummy1_SubField = U_Field.getSubset(
                region=dummy1_nSet, position=NODAL
            )
            U_dummy2_SubField = U_Field.getSubset(
                region=dummy2_nSet, position=NODAL
            )
            #
            if isinstance(U_dummy1_SubField.values[0].data, float):
                # Then variable is a scalar:
                U_dummy1[iFrame] = U_dummy1_SubField.values[0].data
                U_dummy2[iFrame] = U_dummy2_SubField.values[0].data
                #
            else:
                # Variable is an array:
                for j in range(0, len(U_dummy1_SubField.values[0].data)):
                    U_dummy1[iFrame][j] = U_dummy1_SubField.values[0].data[
                        j
                    ]
                    U_dummy2[iFrame][j] = U_dummy2_SubField.values[0].data[
                        j
                    ]

            # Finished saving this variable at the dummy nodes!
            # Variable: RF
            RF_Field = RVEframe.fieldOutputs["RF"]
            RF_dummy1_SubField = RF_Field.getSubset(
                region=dummy1_nSet, position=NODAL
            )
            RF_dummy2_SubField = RF_Field.getSubset(
                region=dummy2_nSet, position=NODAL
            )
            #
            if isinstance(RF_dummy1_SubField.values[0].data, float):
                # Then variable is a scalar:
                RF_dummy1[iFrame] = RF_dummy1_SubField.values[0].data
                RF_dummy2[iFrame] = RF_dummy2_SubField.values[0].data
                #
            else:
                # Variable is an array:
                for j in range(0, len(RF_dummy1_SubField.values[0].data)):
                    RF_dummy1[iFrame][j] = RF_dummy1_SubField.values[
                        0
                    ].data[j]
                    RF_dummy2[iFrame][j] = RF_dummy2_SubField.values[
                        0
                    ].data[j]
            # get the von mises stress
            S_Field = RVEframe.fieldOutputs["S"]
            S_SubField = S_Field.getSubset(
                region=entireRVE_elSet, position=INTEGRATION_POINT
            )

            if isinstance(S_SubField.values[0].mises, float):
                # Then variable is a scalar:
                print("mises  stress is a scalar for each element")
                for i in range(0, len(S_SubField.values)):
                    VonMises[iFrame] = VonMises[iFrame] + \
                        S_SubField.values[i].mises
                # average the value of this increment
                VonMises[iFrame] = VonMises[iFrame] / len(S_SubField.values)
            else:
                # Variable is an array:
                # Loop over every element to compute average
                for j in range(0, len(S_SubField.values[0].mises)):
                    for i in range(0, len(S_SubField.values)):
                        VonMises[iFrame][j] = (
                            VonMises[iFrame][j] + S_SubField.values[i].mises[j]
                        )
                    # average the value of this increment
                    VonMises[iFrame][j] = VonMises[iFrame][j] / len(
                        S_SubField.values
                    )
            # get the equivalent plastic strain
            PEEQ_Field = RVEframe.fieldOutputs["PEEQ"]
            PEEQ_SubField = PEEQ_Field.getSubset(
                region=entireRVE_elSet, position=INTEGRATION_POINT
            )
            if isinstance(PEEQ_SubField.values[0].data, float):
                # Then variable is a scalar:
                print("PEEQ is a scalar for each element")
                for i in range(0, len(PEEQ_SubField.values)):
                    PEEQ[iFrame] = PEEQ[iFrame] + PEEQ_SubField.values[i].data
                # average the value of this increment
                PEEQ[iFrame] = PEEQ[iFrame] / len(PEEQ_SubField.values)
            else:
                # Variable is an array:
                # Loop over every element to compute average
                for j in range(0, len(PEEQ_SubField.values[0].data)):
                    for i in range(0, len(PEEQ_SubField.values)):
                        PEEQ[iFrame][j] = (
                            PEEQ[iFrame][j] + PEEQ_SubField.values[i].data[j]
                        )
                    # average the value of this increment
                    PEEQ[iFrame][j] = PEEQ[iFrame][j] / len(
                        PEEQ_SubField.values)

            # Finished saving this variable at the dummy nodes!
            # Now compute the deformation gradient, jacobian and nominal stress:
            for j in range(0, 2):
                deformation_gradient[iFrame][0][j] = U_dummy1[iFrame][j] + eye[0][j]
                deformation_gradient[iFrame][1][j] = U_dummy2[iFrame][j] + eye[1][j]
                strain[iFrame][0][j] = U_dummy1[iFrame][j]
                strain[iFrame][1][j] = U_dummy2[iFrame][j]
                nominal_stress[iFrame][0][j] = RF_dummy1[iFrame][j] / tot_vol
                nominal_stress[iFrame][1][j] = RF_dummy2[iFrame][j] / tot_vol
    # Save all variables to a single structured variable with all the data
    results = {
        "total_time": total_time,
        "stress": nominal_stress,
        "deformation_gradient": deformation_gradient,
        "strain": strain,
        "total_volume": tot_vol,
        "plastic_energy": ALLPD,
        "VonMises": VonMises,
        "PEEQ": PEEQ,
    }

    # Save post-processing information to pkl file:
    with open("results.pkl", "wb") as fp:
        pickle.dump(results, fp)
