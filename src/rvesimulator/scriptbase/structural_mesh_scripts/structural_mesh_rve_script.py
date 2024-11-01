# -*- coding: mbcs -*-
import numpy as np
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

    # material properties
    youngs_modulus_fiber = sim_info["youngs_modulus_fiber"]
    poisson_ratio_fiber = sim_info["poisson_ratio_fiber"]
    hardening_table_fiber = np.zeros(
        (len(sim_info["hardening_table_fiber"][0]), 2)
    )
    for ii in range(2):
        hardening_table_fiber[:, ii] = sim_info["hardening_table_fiber"][ii]

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
    material_fiber.Plastic(table=(hardening_table_fiber))
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
            'SDV'
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
