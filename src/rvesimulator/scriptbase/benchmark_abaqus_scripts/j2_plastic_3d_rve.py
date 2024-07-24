

try:
    import cPickle as pickle  # Improve speed
except:
    import pickle

import json
from scipy.spatial import distance
from scipy.optimize import linear_sum_assignment

import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
from odbAccess import *

# create a new model
def create_sphere(model, sphere_center, radius, index):
    # center 
    # create a sphere
    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=1.0)
    # create a construction line
    sketch.ConstructionLine(point1=(0.0, -0.5), point2=(0.0, 0.5))
    sketch.FixedConstraint(entity=sketch.geometry[2])
    # information of the sphere (create it at the (0, 0, 0))
    sketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -radius), point2=(0.0, radius), 
        direction=CLOCKWISE)
    # create a line to close the sketch
    sketch.Line(point1=(0.0, -radius), point2=(0.0, radius))
    # revolve the sketch to create a sphere (360 degree)
    sketch.VerticalConstraint(entity=sketch.geometry[4], addUndoState=False)
    # create a new part for the sphere
    part = model.Part(name=f'Sphere-{index}', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    part.BaseSolidRevolve(sketch=sketch, angle=360.0, flipRevolveDirection=OFF)
    # divide the part into four subparts for meshing
    part.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=YZPLANE)
    # create a partition
    part.PartitionCellByDatumPlane(datumPlane=part.datums[2], cells=part.cells[:])
    part.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
    # create a partition
    part.PartitionCellByDatumPlane(datumPlane=part.datums[4], cells=part.cells[:])
    part.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XYPLANE)
    # create a partition
    part.PartitionCellByDatumPlane(datumPlane=part.datums[6], cells=part.cells[:])


    del sketch
    assembly = model.rootAssembly
    instance = assembly.Instance(name=f'Sphere-{index}', part=part, dependent=ON)
    assembly.translate(instanceList=(f'Sphere-{index}', ), vector=sphere_center)

    return instance

def transfer_node_to_array(nodes):
    # transfer the nodes to a numpy array
    nodes_coordinates = np.zeros((len(nodes), 3))
    for ii in range(len(nodes)):
        nodes_coordinates[ii] = np.array(nodes[ii].coordinates)
    return nodes_coordinates

def find_node_pair(face_nodes_list_1, face_nodes_list_2):

    face_nodes_array_1 = transfer_node_to_array(face_nodes_list_1)
    face_nodes_array_2 = transfer_node_to_array(face_nodes_list_2)
    # find the distance between the nodes
    distance_matrix = distance.cdist(face_nodes_array_1, face_nodes_array_2,  'euclidean')
    # find the minimum distance pair
    row_ind, col_ind = linear_sum_assignment(distance_matrix)
    mapping_pairs = list(zip(row_ind, col_ind))
    return mapping_pairs 
def get_node_label(node):
    return node.label

def get_node_x(node):
    return node.coordinates[0]
    
def get_node_y(node):
    return node.coordinates[1]

def get_node_z(node):
    return node.coordinates[2]

def create_random_path(model, path_table, path_name):
    """create a tabular amplitude for history dependent loading path

    Parameters
    ----------
    path_table : list
        a table contains the amplitude
    path_name : str
        name of the amplitude
    """

    # generate the table
    model.TabularAmplitude(
        name=path_name,
        timeSpan=TOTAL,
        smooth=SOLVER_DEFAULT,
        data=(path_table),
    )

def j2_plastic_3d_rve(dict):
    # read the json file using python 2.7 by default

    sphere_info = dict["location_information"]
    len_start = dict["len_start"]
    len_end = dict["len_end"]
    wid_start = dict["wid_start"]
    wid_end = dict["wid_end"]
    height_start = dict["hei_start"]
    height_end = dict["hei_end"]
    radius_mu = dict["radius_mu"]
    radius_std = dict["radius_std"]

    # path dependent information 
    if "strain_amplitude" in dict.keys():
        strain_amplitude = np.zeros(
            (len(dict["strain_amplitude"][0]), 6)
        )
        for ii in range(6):
            strain_amplitude[:, ii] = dict["strain_amplitude"][ii]
        print("path dependent load is applied")
    else:
        print("regular load is applied ")
        strain_amplitude = None

    # material properties (j2 plasticity for both matrix and fibre)
    # matrix
    young_modulus_matrix = dict["youngs_modulus_matrix"]
    poisson_ratio_matrix = dict["poisson_ratio_matrix"]
    hardening_table_matrix = np.zeros((len(dict["hardening_table_matrix"][0]), 2))
    for ii in range(2):
        hardening_table_matrix[:, ii] = dict["hardening_table_matrix"][ii]

    # fiber
    young_modulus_fibre = dict["youngs_modulus_fiber"]
    poisson_ratio_fibre = dict["poisson_ratio_fiber"]
    hardening_table_fibre = np.zeros((len(dict["hardening_table_fiber"][0]), 2))
    for ii in range(2):
        hardening_table_fibre[:, ii] = dict["hardening_table_fiber"][ii]

    # simulation parameters
    time_period = dict["simulation_time"]
    time_interval = dict["simulation_time"]/dict["num_steps"]
    mesh_partition = dict["mesh_partition"]
    
    # loading 
    strain = dict["strain"]
    # information of the RVE
    length = len_end - len_start - 2 * radius_mu
    width = wid_end - wid_start - 2 * radius_mu
    height = height_end - height_start - 2 * radius_mu

    # center of the RVE
    RVEcenter = [0.5*(len_end - len_start - 2*radius_mu), 
                0.5*(wid_end -  wid_start - 2*radius_mu), 
                0.5*(height_end - height_start - 2*radius_mu)]

    delta = 1e-6
    Mesh_size = min(length, width, height)/mesh_partition
    num_cpus = int(dict["num_cpu"])

    # general information 
    model_name = "3d_rve"
    job_name = str(dict["job_name"])

    # create the new model
    mdb.models.changeKey(fromName='Model-1', toName=model_name)
    model = mdb.models[model_name]
    # delete existed model
    if 'Model-1' in mdb.models.keys():
        del mdb.models['Model-1']

    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=1.0)

    sketch.rectangle(point1=(0.0, 0.0), point2=(length, width))
    part = model.Part(name='Matrix', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidExtrude(sketch=sketch, depth=height)
    del sketch
    # create a new assembly for the matrix
    assembly = model.rootAssembly
    instance = assembly.Instance(name='Matrix', part=part, dependent=ON)

    # have a tuple will instance of matrix 
    instance_temp = (instance, )
    # create the spheres and translate them into the corrsponding locations
    for ii  in range(len(sphere_info)):
        sphere_center = (sphere_info[ii][0], sphere_info[ii][1], sphere_info[ii][2])
        shpere_instance = create_sphere(model, sphere_center, sphere_info[ii][3], ii)
        # append the instance to a tuple
        instance_temp = instance_temp + (shpere_instance, )

    # boolean operation to merge the spheres with matrix 
    assembly.InstanceFromBooleanMerge(name="RVE_temp", instances = instance_temp ,
        keepIntersections=ON,  originalInstances=SUPPRESS, domain=GEOMETRY)
    # change name of the instance
    assembly.features.changeKey(fromName='RVE_temp-1',  toName='RVE_temp')

    # create a new part to cut the RVE
    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=1.0)
    sketch.rectangle(point1=(0.0, 0.0), point2=(length, width))
    sketch.rectangle(point1=(-length,-width), point2=(2*length, 2*width))
    part = model.Part(name='Matrix_cut_xoy', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidExtrude(sketch=sketch, depth=2*height)
    del sketch
    # create a new assembly for the matrix
    assembly = model.rootAssembly
    instance = assembly.Instance(name='Matrix_cut_xoy', part=part, dependent=ON)
    # translate the part to the right position 
    assembly.translate(instanceList=(f'Matrix_cut_xoy', ), vector=(0, 0, -0.5*height))
    assembly.InstanceFromBooleanCut(name="RVE_cut", instanceToBeCut=assembly.instances['RVE_temp'],
        cuttingInstances=(assembly.instances['Matrix_cut_xoy'], ), originalInstances=SUPPRESS)
    assembly.features.changeKey(fromName='RVE_cut-1',  toName='RVE_cut')
    # use the cutting part again to cut the RVE
    assembly.Instance(name='Matrix_cut_xoy-rotate', part=part, dependent=ON)
    assembly.rotate(instanceList=('Matrix_cut_xoy-rotate', ), axisPoint=(0.0, 0.0, 0.0), 
        axisDirection=(1.0, 0.0, 0.0), angle=90.0)
    # translate the part to the right position
    assembly.translate(instanceList=('Matrix_cut_xoy-rotate', ), vector=(0.0, 1.5*height, 0.0))

    # cut the RVE again
    assembly.InstanceFromBooleanCut(name="RVE", instanceToBeCut=assembly.instances['RVE_cut'],
        cuttingInstances=(assembly.instances['Matrix_cut_xoy-rotate'], ), originalInstances=SUPPRESS)
    assembly.features.changeKey(fromName='RVE-1',  toName='RVE')

    # delete the unnecessary parts and instances
    for ii, key in enumerate(assembly.features.keys()):
        if key  != 'RVE':
            del assembly.features[key]

    for ii, key in enumerate(model.parts.keys()):
        if key  != 'RVE':
            del model.parts[key]
    part  = model.parts['RVE']

    # create the sets for the RVE
    c = part.cells
    part.Set(cells=c[:], name='All_cells')
    # goes to the cells and find the cell with largest volume
    max_volume = 0
    max_volume_cell_index = 0
    for ii in range(len(c)):
        volume = c[ii].getSize()
        if volume > max_volume:
            max_volume = volume
            max_volume_cell_index = ii
    # create a set for the cell with exculde largest volume
    part.Set(cells=(c[max_volume_cell_index:max_volume_cell_index+1],), name='Matrix')
    # the remaining cells are the spheres using boolean operation
    part.SetByBoolean(
                name="Fibres",
                sets=(part.sets["All_cells"], part.sets["Matrix"]),
                operation=DIFFERENCE,
            )

    # # create sets for every faces
    f = part.faces
    facesBACK = f.getByBoundingBox(RVEcenter[0]-length/2-delta, 
                                RVEcenter[1]-width/2-delta,
                                RVEcenter[2]-height/2-delta, 
                                RVEcenter[0]-length/2+delta, 
                                RVEcenter[1]+width/2+delta, 
                                RVEcenter[2]+height/2+delta)
    part.Set(faces=facesBACK, name='facesBACK')
    facesFRONT = f.getByBoundingBox(RVEcenter[0]+length/2-delta, 
                                    RVEcenter[1]-width/2-delta,
                                    RVEcenter[2]-height/2-delta, 
                                    RVEcenter[0]+length/2+delta, 
                                    RVEcenter[1] + width/2+delta, 
                                    RVEcenter[2]+height/2+delta)
    part.Set(faces=facesFRONT, name='facesFRONT')
    facesTOP = f.getByBoundingBox(RVEcenter[0]-length/2-delta, 
                                RVEcenter[1]-width/2-delta,
                                RVEcenter[2]+height/2-delta, 
                                RVEcenter[0] + length/2+delta, 
                                RVEcenter[1]+width/2+delta,
                                RVEcenter[2]+height/2+delta)
    part.Set(faces=facesTOP, name='facesTOP')
    facesBOTTOM = f.getByBoundingBox(RVEcenter[0]-length/2-delta, 
                                    RVEcenter[1]-width/2-delta,
                                    RVEcenter[2]-height/2-delta, 
                                    RVEcenter[0]+length/2+delta, 
                                    RVEcenter[1]+width/2+delta, 
                                    RVEcenter[2]-height/2+delta)
    part.Set(faces=facesBOTTOM, name='facesBOTTOM')
    facesLEFT = f.getByBoundingBox(RVEcenter[0]-length/2-delta, 
                                RVEcenter[1]-width/2-delta,
                                RVEcenter[2]-height/2-delta, 
                                RVEcenter[0] + length/2+delta, 
                                RVEcenter[1]-width/2+delta, 
                                RVEcenter[2]+height/2+delta)
    part.Set(faces=facesLEFT, name='facesLEFT')
    facesRIGHT = f.getByBoundingBox(RVEcenter[0]-length/2-delta, 
                                    RVEcenter[1]+width/2-delta,
                                    RVEcenter[2]-height/2-delta, 
                                    RVEcenter[0]+length/2+delta,
                                    RVEcenter[1]+width/2+delta, 
                                    RVEcenter[2]+height/2+delta)
    part.Set(faces=facesRIGHT, name='facesRIGHT')

    # create sets for every edges
    s = part.edges
    edgesFRONT_RIGHT = s.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                        RVEcenter[1]+width/2-delta, 
                                        RVEcenter[2]-height/2-delta,
                                            RVEcenter[0]+length/2+delta,
                                            RVEcenter[1]+width/2+delta,
                                                RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesFRONT_RIGHT, name='edgesFRONT_RIGHT')
    edgesFRONT_LEFT = s.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                        RVEcenter[1]-width/2-delta,
                                            RVEcenter[2]-height/2-delta, 
                                            RVEcenter[0]+length/2+delta,
                                            RVEcenter[1]-width/2+delta,
                                            RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesFRONT_LEFT, name='edgesFRONT_LEFT')
    edgesFRONT_TOP = s.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                        RVEcenter[1]-width/2-delta, 
                                        RVEcenter[2]+height/2-delta, 
                                        RVEcenter[0]+length/2+delta, 
                                        RVEcenter[1]+width/2+delta, 
                                        RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesFRONT_TOP, name='edgesFRONT_TOP')
    edgesFRONT_BOTTOM = s.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                            RVEcenter[1]-width/2-delta,
                                            RVEcenter[2]-height/2-delta,
                                            RVEcenter[0]+length/2+delta,
                                            RVEcenter[1]+width/2+delta,
                                                RVEcenter[2]-height/2+delta)
    part.Set(edges=edgesFRONT_BOTTOM, name='edgesFRONT_BOTTOM')
    edgesBACK_RIGHT = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]+width/2-delta, 
                                        RVEcenter[2]-height/2-delta, 
                                        RVEcenter[0]-height/2+delta, 
                                        RVEcenter[1]+width/2+delta,
                                        RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesBACK_RIGHT, name='edgesBACK_RIGHT')
    edgesBACK_LEFT = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]-width/2-delta,
                                        RVEcenter[2]-height/2-delta,
                                        RVEcenter[0]-length/2+delta,
                                        RVEcenter[1]-width/2+delta,
                                        RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesBACK_LEFT, name='edgesBACK_LEFT')
    edgesBACK_TOP = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]-width/2-delta,
                                        RVEcenter[2]+height/2-delta,
                                        RVEcenter[0]-length/2+delta,
                                        RVEcenter[1]+width/2+delta,
                                        RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesBACK_TOP, name='edgesBACK_TOP')
    edgesBACK_BOTTOM = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]-width/2-delta,
                                        RVEcenter[2]-height/2-delta,
                                        RVEcenter[0]-length/2+delta,
                                        RVEcenter[1]+width/2+delta,
                                        RVEcenter[2]-height/2+delta)
    part.Set(edges=edgesBACK_BOTTOM, name='edgesBACK_BOTTOM')
    edgesLEFT_BOTTOM = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]-width/2-delta,
                                        RVEcenter[2]-height/2-delta,
                                        RVEcenter[0]+length/2+delta,
                                        RVEcenter[1]-width/2+delta,
                                        RVEcenter[2]-height/2+delta)
    part.Set(edges=edgesLEFT_BOTTOM, name='edgesLEFT_BOTTOM')
    edgesLEFT_TOP = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]-width/2-delta,
                                        RVEcenter[2]+height/2-delta,
                                        RVEcenter[0]+length/2+delta,
                                        RVEcenter[1]-width/2+delta,
                                        RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesLEFT_TOP, name='edgesLEFT_TOP')
    edgesRIGHT_TOP = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]+width/2-delta,
                                        RVEcenter[2]+height/2-delta,
                                        RVEcenter[0]+length/2+delta, 
                                        RVEcenter[1]+width/2+delta,
                                        RVEcenter[2]+height/2+delta)
    part.Set(edges=edgesRIGHT_TOP, name='edgesRIGHT_TOP')
    edgesRIGHT_BOTTOM = s.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]+width/2-delta,
                                            RVEcenter[2]-height/2-delta, 
                                            RVEcenter[0]+length/2+delta,
                                            RVEcenter[1]+width/2+delta,
                                                RVEcenter[2]-height/2+delta)
    part.Set(edges=edgesRIGHT_BOTTOM, name='edgesRIGHT_BOTTOM')

    # Create sets for every vertices
    v = part.vertices
    vertex_FRB = v.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                    RVEcenter[1]+width/2- delta, 
                                    RVEcenter[2]-height/2-delta, 
                                    RVEcenter[0]+length/2+delta,
                                    RVEcenter[1]+width/2+delta,
                                    RVEcenter[2]-height/2+delta)
    part.Set(vertices=vertex_FRB, name='vertex_FRB')
    vertex_FRT = v.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                    RVEcenter[1]+width/2-delta, 
                                    RVEcenter[2]+height/2-delta, 
                                    RVEcenter[0]+length/2+delta,
                                    RVEcenter[1]+width/2+delta,
                                    RVEcenter[2]+height/2+delta)
    part.Set(vertices=vertex_FRT, name='vertex_FRT')
    vertex_FLT = v.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                    RVEcenter[1]-width/2-delta,
                                    RVEcenter[2]+height/2-delta, 
                                    RVEcenter[0]+length/2+delta,
                                    RVEcenter[1]-width/2+delta,
                                    RVEcenter[2]+height/2+delta)
    part.Set(vertices=vertex_FLT, name='vertex_FLT')
    vertex_FLB = v.getByBoundingBox(RVEcenter[0]+length/2-delta, 
                                    RVEcenter[1]-width/2-delta,
                                    RVEcenter[2]-height/2-delta, 
                                    RVEcenter[0]+length/2+delta,
                                    RVEcenter[1]-width/2+delta,
                                    RVEcenter[2]-height/2+delta)
    part.Set(vertices=vertex_FLB, name='vertex_FLB')
    vertex_BLB = v.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                    RVEcenter[1]-width/2-delta,
                                    RVEcenter[2]-height/2-delta,
                                    RVEcenter[0]-length/2+delta,
                                    RVEcenter[1]-width/2+delta,
                                    RVEcenter[2]-height/2+delta)
    part.Set(vertices=vertex_BLB, name='vertex_BLB')
    vertex_BLT = v.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                    RVEcenter[1]-width/2-delta,
                                    RVEcenter[2]+height/2-delta,
                                    RVEcenter[0]-length/2+delta,
                                    RVEcenter[1]-width/2+delta, 
                                    RVEcenter[2]+height/2+delta)
    part.Set(vertices=vertex_BLT, name='vertex_BLT')
    vertex_BRT = v.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                    RVEcenter[1]+width/2-delta,
                                    RVEcenter[2]+height/2-delta,
                                    RVEcenter[0]-length/2+delta,
                                    RVEcenter[1]+width/2+delta,
                                    RVEcenter[2]+height/2+delta)
    part.Set(vertices=vertex_BRT, name='vertex_BRT')
    vertex_BRB = v.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                    RVEcenter[1]+width/2-delta,
                                    RVEcenter[2]-height/2-delta, 
                                    RVEcenter[0]-length/2+delta,
                                    RVEcenter[1]+width/2+delta,
                                    RVEcenter[2]-height/2+delta)
    part.Set(vertices=vertex_BRB, name='vertex_BRB')

    ## mesh for the RVE ==========================================================
    # set all regions to be tet mesh
    pickedRegions = part.cells[:]
    part.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
    part.seedPart(size=Mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    part.generateMesh()

    # ## define the material and assign material sections
    # 
    model.Material(name='fiber')
    model.materials['fiber'].Elastic(table=((young_modulus_fibre, poisson_ratio_fibre), ))
    model.materials['fiber'].Plastic(table=(hardening_table_fibre))
    model.HomogeneousSolidSection(name='fiber', material='fiber', thickness=None)
    part.SectionAssignment(region=part.sets['Fibres'], sectionName='fiber', offset=0.0,
                        offsetType= MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
    # define the material for the matrix
    model.Material(name='matrix')
    model.materials['matrix'].Elastic(table=((young_modulus_matrix, poisson_ratio_matrix), ))
    model.materials['matrix'].Plastic(table=(hardening_table_matrix))
    model.HomogeneousSolidSection(name='matrix', material='matrix', thickness=None)
    part.SectionAssignment(region=part.sets['Matrix'], sectionName='matrix', offset=0.0,
                        offsetType= MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
    # # pbs: trying to realize the pbc
    assembly.regenerate()
    # define the dummy node by using the reference points
    RF_X_id = assembly.ReferencePoint(point=(RVEcenter[0] + length,RVEcenter[1], RVEcenter[2])).id
    RF_Y_id = assembly.ReferencePoint(point=(RVEcenter[0], RVEcenter[1] + width, RVEcenter[2] )).id
    RF_Z_id = assembly.ReferencePoint(point=(RVEcenter[0], RVEcenter[1], RVEcenter[2]+ height)).id
    refpoints = assembly.referencePoints
    assembly.Set(name='Ref-X', referencePoints=((refpoints[RF_X_id],)))
    assembly.Set(name='Ref-Y', referencePoints=((refpoints[RF_Y_id],)))
    assembly.Set(name='Ref-Z', referencePoints=((refpoints[RF_Z_id],)))

    # find out the Vertices
    vertex_FRB = part.sets['vertex_FRB'].nodes
    assembly.SetFromNodeLabels(name='vertex_FRB', nodeLabels=(('RVE', (vertex_FRB[0].label,)),), unsorted=True)
    vertex_FRT = part.sets['vertex_FRT'].nodes
    assembly.SetFromNodeLabels(name='vertex_FRT', nodeLabels=(('RVE', (vertex_FRT[0].label,)),), unsorted=True)
    vertex_FLT = part.sets['vertex_FLT'].nodes
    assembly.SetFromNodeLabels(name='vertex_FLT', nodeLabels=(('RVE', (vertex_FLT[0].label,)),), unsorted=True)
    vertex_FLB = part.sets['vertex_FLB'].nodes
    assembly.SetFromNodeLabels(name='vertex_FLB', nodeLabels=(('RVE', (vertex_FLB[0].label,)),), unsorted=True)
    vertex_BLB = part.sets['vertex_BLB'].nodes
    assembly.SetFromNodeLabels(name='vertex_BLB', nodeLabels=(('RVE', (vertex_BLB[0].label,)),), unsorted=True)
    vertex_BLT = part.sets['vertex_BLT'].nodes
    assembly.SetFromNodeLabels(name='vertex_BLT', nodeLabels=(('RVE', (vertex_BLT[0].label,)),), unsorted=True)
    vertex_BRT = part.sets['vertex_BRT'].nodes
    assembly.SetFromNodeLabels(name='vertex_BRT', nodeLabels=(('RVE', (vertex_BRT[0].label,)),), unsorted=True)
    vertex_BRB = part.sets['vertex_BRB'].nodes
    assembly.SetFromNodeLabels(name='vertex_BRB', nodeLabels=(('RVE', (vertex_BRB[0].label,)),), unsorted=True)

    # define the pbc for vertices ========================================
    # part 1: equations for vertices 3(FRT) and 5(BLB)
    for ii in range(1, 4):
        model.Equation(name='FRT_BLB_'+str(ii), 
                    terms=((1, 'vertex_FRT', ii),
                            (-1, 'vertex_BLB', ii),
                            (-1*length, 'Ref-X', ii),
                            (-1*width, 'Ref-Y', ii),
                            (-1*height, 'Ref-Z', ii)))
    # part 2: equations for vertices 2(FRB) and 8 (BLT)
    for ii in range(1, 4):
        model.Equation(name='FRB_BLT_'+str(ii), 
                    terms=((1, 'vertex_FRB', ii), 
                            (-1, 'vertex_BLT', ii), 
                            (-1*length, 'Ref-X', ii), 
                            (-1*width, 'Ref-Y', ii), 
                            (1*height, 'Ref-Z', ii)))
    # part 3: equations for vertices 7 (BRT) and 1 (FLB)
    for ii in range(1, 4):
        model.Equation(name='BRT_FLB_'+str(ii),
                        terms=((1, 'vertex_BRT', ii),
                            (-1, 'vertex_FLB', ii),
                            (1*length, 'Ref-X', ii), 
                            (-1*width, 'Ref-Y', ii),
                            (-1*height, 'Ref-Z', ii)))
    # part 4: equations for vertices 4 (BRT) and 6 (FLB)
    for ii in range(1, 4):
        model.Equation(name='FLT_BRB_'+str(ii), 
                    terms=((1, 'vertex_FLT', ii),
                            (-1, 'vertex_BRB', ii),
                            (-1*length, 'Ref-X', ii),
                            (1*width, 'Ref-Y', ii),
                            (-1*height, 'Ref-Z', ii)))

    ## define the pbc for edges ========================================

    # part 1: equations for edges 2 (edgesFRONT_RIGHT) and 4 (edgesBACK_LEFT)
    edgesFRONT_RIGHT_nodes = part.sets['edgesFRONT_RIGHT'].nodes
    edgesFRONT_RIGHT_nodes_sorted = sorted(edgesFRONT_RIGHT_nodes, key=get_node_z)
    # pop the first and last node (those are the vertices) 
    edgesFRONT_RIGHT_nodes_sorted.pop(0)
    edgesFRONT_RIGHT_nodes_sorted.pop(-1)
    edgesBACK_LEFT_nodes = part.sets['edgesBACK_LEFT'].nodes
    edgesBACK_LEFT_nodes_sorted = sorted(edgesBACK_LEFT_nodes, key=get_node_z)
    # pop the first and last node (those are the vertices)
    edgesBACK_LEFT_nodes_sorted.pop(0)
    edgesBACK_LEFT_nodes_sorted.pop(-1)
    # check if the number of nodes are the same
    if len(edgesFRONT_RIGHT_nodes_sorted) ==len(edgesBACK_LEFT_nodes_sorted):
        for ii in range(len(edgesFRONT_RIGHT_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTRIGHT_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_RIGHT_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKLEFT_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_LEFT_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1,4):
                model.Equation(name='FRONTRIGHT_BACKLEFT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTRIGHT_' + str(ii), jj), (-1, 'BACKLEFT_' + str(ii), jj), (-1*length, 'Ref-X', jj),(-1*width, 'Ref-Y', jj)))
    else:
        print ("The number of nodes between the two sides are not the same")
        print (f"edgesFRONT_RIGHT_nodes: {len(edgesFRONT_RIGHT_nodes_sorted)}, edgesBACK_LEFT_nodes:{len(edgesBACK_LEFT_nodes_sorted)}")
        print (f"Randomly to remove {np.abs(len(edgesFRONT_RIGHT_nodes_sorted) - len(edgesBACK_LEFT_nodes_sorted))} nodes from the larger side")
        if len(edgesFRONT_RIGHT_nodes_sorted) > len(edgesBACK_LEFT_nodes_sorted):
            num_remove = len(edgesFRONT_RIGHT_nodes_sorted) - len(edgesBACK_LEFT_nodes_sorted)
            # randomly generate the index to remove
            index_remove = np.random.choice(len(edgesFRONT_RIGHT_nodes_sorted)- num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesFRONT_RIGHT_nodes_sorted.pop(index_remove[jj])
        
        else:
            num_remove = len(edgesBACK_LEFT_nodes_sorted) - len(edgesFRONT_RIGHT_nodes_sorted)
            index_remove = np.random.choice(len(edgesBACK_LEFT_nodes_sorted)- num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesBACK_LEFT_nodes_sorted.pop(index_remove[jj])
        
        # apply the pbc for the nodes after removing the nodes
        for ii in range(len(edgesFRONT_RIGHT_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTRIGHT_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_RIGHT_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKLEFT_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_LEFT_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1,4):
                model.Equation(name='FRONTRIGHT_BACKLEFT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTRIGHT_' + str(ii), jj), (-1, 'BACKLEFT_' + str(ii), jj), (-1*length, 'Ref-X', jj),(-1*width, 'Ref-Y', jj)))
            

    # part 2: equations for edges 1 (edgesFRONT_LEFT) and 3 (edgesBACK_RIGHT)
    edgesFRONT_LEFT_nodes = part.sets['edgesFRONT_LEFT'].nodes
    edgesFRONT_LEFT_nodes_sorted = sorted(edgesFRONT_LEFT_nodes, key=get_node_z)
    edgesFRONT_LEFT_nodes_sorted.pop(0)
    edgesFRONT_LEFT_nodes_sorted.pop(-1)

    edgesBACK_RIGHT_nodes = part.sets['edgesBACK_RIGHT'].nodes
    edgesBACK_RIGHT_nodes_sorted = sorted(edgesBACK_RIGHT_nodes, key=get_node_z)
    edgesBACK_RIGHT_nodes_sorted.pop(0)
    edgesBACK_RIGHT_nodes_sorted.pop(-1)

    if len(edgesFRONT_LEFT_nodes_sorted) == len(edgesBACK_RIGHT_nodes_sorted):
        for ii in range(len(edgesFRONT_LEFT_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTLEFT_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_LEFT_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKRIGHT_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_RIGHT_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONTLEFT_BACKRIGHT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTLEFT_' + str(ii), jj), (-1, 'BACKRIGHT_' + str(ii), jj), (-1*length, 'Ref-X', jj), (1*width, 'Ref-Y', jj)))
    else:
        print ("The number of nodes between the two sides are not the same")
        print (f"edgesFRONT_LEFT_nodes: {len(edgesFRONT_LEFT_nodes_sorted)}, edgesBACK_RIGHT_nodes:{len(edgesBACK_RIGHT_nodes_sorted)}")
        print (f"Randomly to remove {np.abs(len(edgesFRONT_LEFT_nodes_sorted) - len(edgesBACK_RIGHT_nodes_sorted))} nodes from the larger side")
        if len(edgesFRONT_LEFT_nodes_sorted) > len(edgesBACK_RIGHT_nodes_sorted):
            num_remove = len(edgesFRONT_LEFT_nodes_sorted) - len(edgesBACK_RIGHT_nodes_sorted)
            # randomly generate the index to remove
            index_remove = np.random.choice(len(edgesFRONT_LEFT_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesFRONT_LEFT_nodes_sorted.pop(index_remove[jj])
        else:
            num_remove = len(edgesBACK_RIGHT_nodes_sorted) - len(edgesFRONT_LEFT_nodes_sorted)
            index_remove = np.random.choice(len(edgesBACK_RIGHT_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesBACK_RIGHT_nodes_sorted.pop(index_remove[jj])

        # apply the pbc for the nodes after removing the nodes
        for ii in range(len(edgesFRONT_LEFT_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTLEFT_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_LEFT_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKRIGHT_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_RIGHT_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONTLEFT_BACKRIGHT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTLEFT_' + str(ii), jj), (-1, 'BACKRIGHT_' + str(ii), jj), (-1*length, 'Ref-X', jj), (1*width, 'Ref-Y', jj)))

    # part 3: equations for edges 6 (edgesFRONT_TOP) and 8 (edgesBACK_BOTTOM)
    edgesFRONT_TOP_nodes = part.sets['edgesFRONT_TOP'].nodes
    edgesFRONT_TOP_nodes_sorted = sorted(edgesFRONT_TOP_nodes, key=get_node_y)
    # pop the first and last node (those are the vertices)
    edgesFRONT_TOP_nodes_sorted.pop(0)
    edgesFRONT_TOP_nodes_sorted.pop(-1)
    edgesBACK_BOTTOM_nodes = part.sets['edgesBACK_BOTTOM'].nodes
    edgesBACK_BOTTOM_nodes_sorted = sorted(edgesBACK_BOTTOM_nodes, key=get_node_y)
    # pop the first and last node (those are the vertices)
    edgesBACK_BOTTOM_nodes_sorted.pop(0)
    edgesBACK_BOTTOM_nodes_sorted.pop(-1)

    if len(edgesFRONT_TOP_nodes_sorted) == len(edgesBACK_BOTTOM_nodes_sorted):
        for ii in range(len(edgesFRONT_TOP_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_TOP_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONTTOP_BACKBOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTTOP_' + str(ii), jj), (-1, 'BACKBOTTOM_' + str(ii), jj), (-1*length, 'Ref-X', jj), (-1*height, 'Ref-Z', jj)))
    else:
        print ("The number of nodes between the two sides are not the same")
        print (f"edgesFRONT_TOP_nodes: {len(edgesFRONT_TOP_nodes_sorted)}, edgesBACK_BOTTOM_nodes:{len(edgesBACK_BOTTOM_nodes_sorted)}")
        print (f"Randomly to remove {np.abs(len(edgesFRONT_TOP_nodes_sorted) - len(edgesBACK_BOTTOM_nodes_sorted))} nodes from the larger side")
        if len(edgesFRONT_TOP_nodes_sorted) > len(edgesBACK_BOTTOM_nodes_sorted):
            num_remove = len(edgesFRONT_TOP_nodes_sorted) - len(edgesBACK_BOTTOM_nodes_sorted)
            # randomly generate the index to remove
            index_remove = np.random.choice(len(edgesFRONT_TOP_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesFRONT_TOP_nodes_sorted.pop(index_remove[jj])
        else:
            num_remove = len(edgesBACK_BOTTOM_nodes_sorted) - len(edgesFRONT_TOP_nodes_sorted)
            index_remove = np.random.choice(len(edgesBACK_BOTTOM_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesBACK_BOTTOM_nodes_sorted.pop(index_remove[jj])
        
        # apply the pbc for the nodes after removing the nodes
        for ii in range(len(edgesFRONT_TOP_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_TOP_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONTTOP_BACKBOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTTOP_' + str(ii), jj), (-1, 'BACKBOTTOM_' + str(ii), jj), (-1*length, 'Ref-X', jj), (-1*height, 'Ref-Z', jj)))
        
    # part 4: equations for edges 5 (edgesFRONT_BOTTOM) and 7 (edgesBACK_TOP)
    edgesFRONT_BOTTOM_nodes = part.sets['edgesFRONT_BOTTOM'].nodes
    edgesFRONT_BOTTOM_nodes_sorted = sorted(edgesFRONT_BOTTOM_nodes, key=get_node_y)
    edgesFRONT_BOTTOM_nodes_sorted.pop(0)
    edgesFRONT_BOTTOM_nodes_sorted.pop(-1)

    edgesBACK_TOP_nodes = part.sets['edgesBACK_TOP'].nodes
    edgesBACK_TOP_nodes_sorted = sorted(edgesBACK_TOP_nodes, key=get_node_y)
    edgesBACK_TOP_nodes_sorted.pop(0)
    edgesBACK_TOP_nodes_sorted.pop(-1)
    # check if the number of nodes are the same
    if len(edgesFRONT_BOTTOM_nodes_sorted) == len(edgesBACK_TOP_nodes_sorted):
        for ii in range(1, len(edgesFRONT_BOTTOM_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_TOP_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONTBOTTOM_BACKTOP_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTBOTTOM_' + str(ii), jj), (-1, 'BACKTOP_' + str(ii), jj), (-1*length, 'Ref-X', jj), (1*height, 'Ref-Z', jj)))
    else:
        print ("The number of nodes between the two sides are not the same")
        print (f"edgesFRONT_BOTTOM_nodes: {len(edgesFRONT_BOTTOM_nodes_sorted)}, edgesBACK_TOP_nodes:{len(edgesBACK_TOP_nodes_sorted)}")
        print (f"Randomly to remove {np.abs(len(edgesFRONT_BOTTOM_nodes_sorted) - len(edgesBACK_TOP_nodes_sorted))} nodes from the larger side")
        if len(edgesFRONT_BOTTOM_nodes_sorted) > len(edgesBACK_TOP_nodes_sorted):
            num_remove = len(edgesFRONT_BOTTOM_nodes_sorted) - len(edgesBACK_TOP_nodes_sorted)
            # randomly generate the index to remove
            index_remove = np.random.choice(len(edgesFRONT_BOTTOM_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesFRONT_BOTTOM_nodes_sorted.pop(index_remove[jj])
        else:
            num_remove = len(edgesBACK_TOP_nodes_sorted) - len(edgesFRONT_BOTTOM_nodes_sorted)
            index_remove = np.random.choice(len(edgesBACK_TOP_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesBACK_TOP_nodes_sorted.pop(index_remove[jj])

        # apply the pbc for the nodes after removing the nodes
        for ii in range(1, len(edgesFRONT_BOTTOM_nodes_sorted)):
            assembly.SetFromNodeLabels(name='FRONTBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesFRONT_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='BACKTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesBACK_TOP_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONTBOTTOM_BACKTOP_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONTBOTTOM_' + str(ii), jj), (-1, 'BACKTOP_' + str(ii), jj), (-1*length, 'Ref-X', jj), (1*height, 'Ref-Z', jj)))


    # part 5: equations for edges 11 (edgesRIGHT_TOP) and 9 (edgesLEFT_BOTTOM)
    edgesRIGHT_TOP_nodes = part.sets['edgesRIGHT_TOP'].nodes
    edgesRIGHT_TOP_nodes_sorted = sorted(edgesRIGHT_TOP_nodes, key=get_node_x)
    edgesRIGHT_TOP_nodes_sorted.pop(0)
    edgesRIGHT_TOP_nodes_sorted.pop(-1)
    edgesLEFT_BOTTOM_nodes = part.sets['edgesLEFT_BOTTOM'].nodes
    edgesLEFT_BOTTOM_nodes_sorted = sorted(edgesLEFT_BOTTOM_nodes, key=get_node_x)
    edgesLEFT_BOTTOM_nodes_sorted.pop(0)
    edgesLEFT_BOTTOM_nodes_sorted.pop(-1)
    # check if the number of nodes are the same
    if len(edgesRIGHT_TOP_nodes_sorted) == len(edgesLEFT_BOTTOM_nodes_sorted):
        for ii in range(len(edgesRIGHT_TOP_nodes_sorted)):
            assembly.SetFromNodeLabels(name='RIGHTTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesRIGHT_TOP_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='LEFTBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesLEFT_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHTTOP_LEFTBOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHTTOP_' + str(ii), jj), (-1, 'LEFTBOTTOM_' + str(ii), jj), (-1*width, 'Ref-Y', jj), (-1*height, 'Ref-Z', jj)))
    else:
        print ("The number of nodes between the two sides are not the same")
        print (f"edgesRIGHT_TOP_nodes: {len(edgesRIGHT_TOP_nodes_sorted)}, edgesLEFT_BOTTOM_nodes:{len(edgesLEFT_BOTTOM_nodes_sorted)}")
        print (f"Randomly to remove {np.abs(len(edgesRIGHT_TOP_nodes_sorted) - len(edgesLEFT_BOTTOM_nodes_sorted))} nodes from the larger side")
        if len(edgesRIGHT_TOP_nodes_sorted) > len(edgesLEFT_BOTTOM_nodes_sorted):
            num_remove = len(edgesRIGHT_TOP_nodes_sorted) - len(edgesLEFT_BOTTOM_nodes_sorted)
            # randomly generate the index to remove
            index_remove = np.random.choice(len(edgesRIGHT_TOP_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesRIGHT_TOP_nodes_sorted.pop(index_remove[jj])
        else:
            num_remove = len(edgesLEFT_BOTTOM_nodes_sorted) - len(edgesRIGHT_TOP_nodes_sorted)
            index_remove = np.random.choice(len(edgesLEFT_BOTTOM_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesLEFT_BOTTOM_nodes_sorted.pop(index_remove[jj])
        # apply the pbc for the nodes after removing the nodes
        for ii in range(len(edgesRIGHT_TOP_nodes_sorted)):
            assembly.SetFromNodeLabels(name='RIGHTTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesRIGHT_TOP_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='LEFTBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesLEFT_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHTTOP_LEFTBOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHTTOP_' + str(ii), jj), (-1, 'LEFTBOTTOM_' + str(ii), jj), (-1*width, 'Ref-Y', jj), (-1*height, 'Ref-Z', jj)))

    # part 6: equations for edges 10 (edgesRIGHT_BOTTOM) and 12 (edgesLEFT_TOP)
    edgesRIGHT_BOTTOM_nodes = part.sets['edgesRIGHT_BOTTOM'].nodes
    edgesRIGHT_BOTTOM_nodes_sorted = sorted(edgesRIGHT_BOTTOM_nodes, key=get_node_x)
    edgesRIGHT_BOTTOM_nodes_sorted.pop(0)
    edgesRIGHT_BOTTOM_nodes_sorted.pop(-1)
    edgesLEFT_TOP_nodes = part.sets['edgesLEFT_TOP'].nodes
    edgesLEFT_TOP_nodes_sorted = sorted(edgesLEFT_TOP_nodes, key=get_node_x)
    edgesLEFT_TOP_nodes_sorted.pop(0)
    edgesLEFT_TOP_nodes_sorted.pop(-1)
    if len(edgesRIGHT_BOTTOM_nodes_sorted) == len(edgesLEFT_TOP_nodes_sorted):
        for ii in range(len(edgesRIGHT_BOTTOM_nodes_sorted)):
            assembly.SetFromNodeLabels(name='RIGHTBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesRIGHT_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='LEFTTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesLEFT_TOP_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHTBOTTOM_LEFTTOP_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHTBOTTOM_' + str(ii), jj), (-1, 'LEFTTOP_' + str(ii), jj), (-1*width, 'Ref-Y', jj), (1*height, 'Ref-Z', jj)))
    else:
        print ("The number of nodes between the two sides are not the same")
        print (f"edgesRIGHT_BOTTOM_nodes: {len(edgesRIGHT_BOTTOM_nodes_sorted)}, edgesLEFT_TOP_nodes:{len(edgesLEFT_TOP_nodes_sorted)}")
        print (f"Randomly to remove {np.abs(len(edgesRIGHT_BOTTOM_nodes_sorted) - len(edgesLEFT_TOP_nodes_sorted))} nodes from the larger side")
        if len(edgesRIGHT_BOTTOM_nodes_sorted) > len(edgesLEFT_TOP_nodes_sorted):
            num_remove = len(edgesRIGHT_BOTTOM_nodes_sorted) - len(edgesLEFT_TOP_nodes_sorted)
            # randomly generate the index to remove
            index_remove = np.random.choice(len(edgesRIGHT_BOTTOM_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesRIGHT_BOTTOM_nodes_sorted.pop(index_remove[jj])
        else:
            num_remove = len(edgesLEFT_TOP_nodes_sorted) - len(edgesRIGHT_BOTTOM_nodes_sorted)
            index_remove = np.random.choice(len(edgesLEFT_TOP_nodes_sorted) - num_remove, num_remove, replace=False)
            for jj in range(len(index_remove)):
                edgesLEFT_TOP_nodes_sorted.pop(index_remove[jj])
        # apply the pbc for the nodes after removing the nodes
        for ii in range(len(edgesRIGHT_BOTTOM_nodes_sorted)):
            assembly.SetFromNodeLabels(name='RIGHTBOTTOM_' + str(ii), nodeLabels=(('RVE', tuple([edgesRIGHT_BOTTOM_nodes_sorted[ii].label])),), unsorted=True)
            assembly.SetFromNodeLabels(name='LEFTTOP_' + str(ii), nodeLabels=(('RVE', tuple([edgesLEFT_TOP_nodes_sorted[ii].label])),), unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHTBOTTOM_LEFTTOP_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHTBOTTOM_' + str(ii), jj), (-1, 'LEFTTOP_' + str(ii), jj), (-1*width, 'Ref-Y', jj), (1*height, 'Ref-Z', jj)))

    ###
    # create sets for every faces
    allnodes = part.nodes
    # node on the back face 
    BACK_nodes = allnodes.getByBoundingBox(RVEcenter[0]-length/2-delta,
                                        RVEcenter[1]-width/2+delta,
                                        RVEcenter[2]-height/2+delta, 
                                        RVEcenter[0]-length/2+delta,
                                            RVEcenter[1]+width/2-delta,
                                            RVEcenter[2]+height/2-delta)
    part.Set(nodes=BACK_nodes, name='BACK_nodes')
    # node on the front face
    FRONT_nodes = allnodes.getByBoundingBox(RVEcenter[0]+length/2-delta,
                                            RVEcenter[1]-width/2+delta,
                                            RVEcenter[2]-height/2+delta,
                                            RVEcenter[0]+length/2+delta,
                                            RVEcenter[1]+width/2-delta,
                                            RVEcenter[2]+height/2-delta)
    part.Set(nodes=FRONT_nodes, name='FRONT_nodes')


    TOP_nodes = allnodes.getByBoundingBox(RVEcenter[0]-length/2+delta,
                                        RVEcenter[1]-width/2+delta, 
                                        RVEcenter[2]+height/2-delta,
                                        RVEcenter[0]+length/2-delta,
                                        RVEcenter[1]+width/2-delta,
                                        RVEcenter[2]+height/2+delta)
    part.Set(nodes=TOP_nodes, name='TOP_nodes')
    BOTTOM_nodes =allnodes.getByBoundingBox(RVEcenter[0]-length/2+delta,
                                            RVEcenter[1]-width/2+delta,
                                            RVEcenter[2]-height/2-delta,
                                            RVEcenter[0]+length/2-delta,
                                            RVEcenter[1]+width/2-delta,
                                            RVEcenter[2]-height/2+delta)
    part.Set(nodes=BOTTOM_nodes, name='BOTTOM_nodes')
    LEFT_nodes = allnodes.getByBoundingBox(RVEcenter[0]-length/2+delta,
                                            RVEcenter[1]-width/2-delta,
                                            RVEcenter[2]-height/2+delta,
                                            RVEcenter[0]+length/2-delta,
                                            RVEcenter[1]-width/2+delta, 
                                            RVEcenter[2]+height/2-delta)
    part.Set(nodes=LEFT_nodes, name='LEFT_nodes')
    RIGHT_nodes = allnodes.getByBoundingBox(RVEcenter[0]-length/2+delta,
                                            RVEcenter[1]+width/2-delta,
                                            RVEcenter[2]-height/2+delta,
                                            RVEcenter[0]+length/2-delta, 
                                            RVEcenter[1]+width/2+delta, 
                                            RVEcenter[2]+height/2-delta)
    part.Set(nodes=RIGHT_nodes, name='RIGHT_nodes')

    # create sets for support nodes
    support = allnodes.getByBoundingBox(RVEcenter[0]-5*Mesh_size,
                                        RVEcenter[1]-5*Mesh_size,
                                        RVEcenter[2]-5*Mesh_size,
                                        RVEcenter[0]+5*Mesh_size, 
                                        RVEcenter[1]+5*Mesh_size, 
                                        RVEcenter[2]+5*Mesh_size)
    part.Set(nodes=support, name='support_nodes')
    support = sorted(support, key=get_node_label)  

    # sort the nodes based on the label 
    BACK_nodes = sorted(BACK_nodes, key=get_node_label)
    FRONT_nodes = sorted(FRONT_nodes, key=get_node_label)

    if len(FRONT_nodes) == len(BACK_nodes):
        # fine the mapping between the two faces 
        mapping_pairs  = find_node_pair(FRONT_nodes, BACK_nodes)
        # apply pbc 
        for ii in range(len(FRONT_nodes)):
            assembly.SetFromNodeLabels(name='FRONT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([FRONT_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='BACK_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([BACK_nodes[mapping_pairs[ii][1]].label])),),
                                    unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONT_BACK_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONT_' + str(ii), jj),
                                        (-1, 'BACK_' + str(ii), jj),
                                        (-1*length, 'Ref-X', jj)))
    elif len(FRONT_nodes) > len(BACK_nodes):
        print ("The number of nodes between the two sides are not the same")
        print (f"FRONT_nodes: {len(FRONT_nodes)}, BACK_nodes:{len(BACK_nodes)}")
        print (f"Randomly to remove {np.abs(len(FRONT_nodes) - len(BACK_nodes))} nodes from the larger side")
        num_remove = len(FRONT_nodes) - len(BACK_nodes)
        # randomly generate the index to remove
        index_remove = np.random.choice(len(FRONT_nodes) - num_remove, num_remove, replace=False)
        for jj in range(len(index_remove)):
            FRONT_nodes.pop(index_remove[jj])
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(FRONT_nodes, BACK_nodes)
        # apply pbc
        for ii in range(len(FRONT_nodes)):
            assembly.SetFromNodeLabels(name='FRONT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([FRONT_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='BACK_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([BACK_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONT_BACK_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONT_' + str(ii), jj),
                                        (-1, 'BACK_' + str(ii), jj),
                                        (-1*length, 'Ref-X', jj)))
    elif len(FRONT_nodes) < len(BACK_nodes):
        print ("The number of nodes between the two sides are not the same")
        print (f"FRONT_nodes: {len(FRONT_nodes)}, BACK_nodes:{len(BACK_nodes)}")
        print (f"Randomly to remove {np.abs(len(FRONT_nodes) - len(BACK_nodes))} nodes from the larger side")
        num_remove = len(BACK_nodes) - len(FRONT_nodes)
        # randomly generate the index to remove
        index_remove = np.random.choice(len(BACK_nodes) - num_remove, num_remove, replace=False)
        for jj in range(len(index_remove)):
            BACK_nodes.pop(index_remove[jj])
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(FRONT_nodes, BACK_nodes)
        # apply pbc
        for ii in range(len(FRONT_nodes)):
            assembly.SetFromNodeLabels(name='FRONT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([FRONT_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='BACK_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([BACK_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='FRONT_BACK_' + str(ii) + '_' + str(jj),
                                terms=((1, 'FRONT_' + str(ii), jj),
                                        (-1, 'BACK_' + str(ii), jj),
                                       (-1*length, 'Ref-X', jj)))
    else: 
        raise ValueError("Something wrong with the number of nodes on the two faces")

    # sort the nodes based on the label
    RIGHT_nodes = sorted(RIGHT_nodes, key=get_node_label)
    LEFT_nodes = sorted(LEFT_nodes, key=get_node_label)
    #  for the RIGHT and LEFT faces
    if len(RIGHT_nodes) == len(LEFT_nodes):
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(RIGHT_nodes, LEFT_nodes)
        # apply pbc
        for ii in range(len(RIGHT_nodes)):
            assembly.SetFromNodeLabels(name='RIGHT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([RIGHT_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='LEFT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([LEFT_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHT_LEFT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHT_' + str(ii), jj),
                                        (-1, 'LEFT_' + str(ii), jj),
                                        (-1*width, 'Ref-Y', jj)))
    elif len(RIGHT_nodes) > len(LEFT_nodes):
        print ("The number of nodes between the two sides are not the same")
        print (f"RIGHT_nodes: {len(RIGHT_nodes)}, LEFT_nodes:{len(LEFT_nodes)}")
        print (f"Randomly to remove {np.abs(len(RIGHT_nodes) - len(LEFT_nodes))} nodes from the larger side")
        num_remove = len(RIGHT_nodes) - len(LEFT_nodes)
        # randomly generate the index to remove
        index_remove = np.random.choice(len(RIGHT_nodes) - num_remove, num_remove, replace=False)
        for jj in range(len(index_remove)):
            RIGHT_nodes.pop(index_remove[jj])
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(RIGHT_nodes, LEFT_nodes)
        # apply pbc
        for ii in range(len(RIGHT_nodes)):
            assembly.SetFromNodeLabels(name='RIGHT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([RIGHT_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='LEFT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([LEFT_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHT_LEFT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHT_' + str(ii), jj),
                                        (-1, 'LEFT_' + str(ii), jj),
                                        (-1*width, 'Ref-Y', jj)))
    elif len(RIGHT_nodes) < len(LEFT_nodes):
        print ("The number of nodes between the two sides are not the same")
        print (f"RIGHT_nodes: {len(RIGHT_nodes)}, LEFT_nodes:{len(LEFT_nodes)}")
        print (f"Randomly to remove {np.abs(len(RIGHT_nodes) - len(LEFT_nodes))} nodes from the larger side")
        num_remove = len(LEFT_nodes) - len(RIGHT_nodes)
        # randomly generate the index to remove
        index_remove = np.random.choice(len(LEFT_nodes) - num_remove, num_remove, replace=False)
        for jj in range(len(index_remove)):
            LEFT_nodes.pop(index_remove[jj])
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(RIGHT_nodes, LEFT_nodes)
        # apply pbc
        for ii in range(len(RIGHT_nodes)):
            assembly.SetFromNodeLabels(name='RIGHT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([RIGHT_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='LEFT_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([LEFT_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='RIGHT_LEFT_' + str(ii) + '_' + str(jj),
                                terms=((1, 'RIGHT_' + str(ii), jj),
                                        (-1, 'LEFT_' + str(ii), jj),
                                            (-1*width, 'Ref-Y', jj)))
    else:
        raise ValueError("Something wrong with the number of nodes on the two faces")

    # sort the nodes based on the label
    TOP_nodes = sorted(TOP_nodes, key=get_node_label)
    BOTTOM_nodes = sorted(BOTTOM_nodes, key=get_node_label)
    # for the TOP and BOTTOM faces
    if len(TOP_nodes) == len(BOTTOM_nodes):
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(TOP_nodes, BOTTOM_nodes)
        # apply pbc
        for ii in range(len(TOP_nodes)):
            assembly.SetFromNodeLabels(name='TOP_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([TOP_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='BOTTOM_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([BOTTOM_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='TOP_BOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'TOP_' + str(ii), jj),
                                        (-1, 'BOTTOM_' + str(ii), jj),
                                        (-1*height, 'Ref-Z', jj)))
    elif len(TOP_nodes) > len(BOTTOM_nodes):
        print ("The number of nodes between the two sides are not the same")
        print (f"TOP_nodes: {len(TOP_nodes)}, BOTTOM_nodes:{len(BOTTOM_nodes)}")
        print (f"Randomly to remove {np.abs(len(TOP_nodes) - len(BOTTOM_nodes))} nodes from the larger side")
        num_remove = len(TOP_nodes) - len(BOTTOM_nodes)
        # randomly generate the index to remove
        index_remove = np.random.choice(len(TOP_nodes) - num_remove, num_remove, replace=False)
        for jj in range(len(index_remove)):
            TOP_nodes.pop(index_remove[jj])
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(TOP_nodes, BOTTOM_nodes)
        # apply pbc
        for ii in range(len(TOP_nodes)):
            assembly.SetFromNodeLabels(name='TOP_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([TOP_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='BOTTOM_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([BOTTOM_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='TOP_BOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'TOP_' + str(ii), jj),
                                        (-1, 'BOTTOM_' + str(ii), jj),
                                         (-1*height, 'Ref-Z', jj)))
    elif len(TOP_nodes) < len(BOTTOM_nodes):
        print ("The number of nodes between the two sides are not the same")
        print (f"TOP_nodes: {len(TOP_nodes)}, BOTTOM_nodes:{len(BOTTOM_nodes)}")
        print (f"Randomly to remove {np.abs(len(TOP_nodes) - len(BOTTOM_nodes))} nodes from the larger side")
        num_remove = len(BOTTOM_nodes) - len(TOP_nodes)
        # randomly generate the index to remove
        index_remove = np.random.choice(len(BOTTOM_nodes) - num_remove, num_remove, replace=False)
        for jj in range(len(index_remove)):
            BOTTOM_nodes.pop(index_remove[jj])
        # find the mapping between the two faces 
        mapping_pairs = find_node_pair(TOP_nodes, BOTTOM_nodes)
        # apply pbc
        for ii in range(len(TOP_nodes)):
            assembly.SetFromNodeLabels(name='TOP_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([TOP_nodes[mapping_pairs[ii][0]].label])),), 
                                    unsorted=True)
            assembly.SetFromNodeLabels(name='BOTTOM_' + str(ii), 
                                    nodeLabels=(('RVE', tuple([BOTTOM_nodes[mapping_pairs[ii][1]].label])),),
                                        unsorted=True)
            for jj in range(1, 4):
                model.Equation(name='TOP_BOTTOM_' + str(ii) + '_' + str(jj),
                                terms=((1, 'TOP_' + str(ii), jj),
                                        (-1, 'BOTTOM_' + str(ii), jj),
                                        (-1*height, 'Ref-Z', jj)))

    else:
        raise ValueError("Something wrong with the number of nodes on the two faces")

    # create step 
    model.StaticStep(
                     initialInc=0.01,
                        maxInc=1.0,
                        maxNumInc=1000,
                        minInc=1e-5,
                        name="Step-1",
                        previous="Initial",
                        timePeriod=time_period,
                        nlgeom=ON,)
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
    ### strain for the first dummy point(the front one E11,E12,E13) 

    if strain_amplitude is None:
    
        model.DisplacementBC(name='E_11', createStepName='Step-1',
                                region=assembly.sets['Ref-X'], u1=strain[0], u2=UNSET, u3=UNSET, amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_12', createStepName='Step-1',
                                region=assembly.sets['Ref-X'], u1=UNSET, u2=strain[3], u3=UNSET, amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_13', createStepName='Step-1',
                                region=assembly.sets['Ref-X'], u1=UNSET, u2=UNSET, u3=strain[5], amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        
        #### strain for the second dummy point(the front one E21,E22,E23)
        model.DisplacementBC(name='E_21', createStepName='Step-1',
                                region=assembly.sets['Ref-Y'], u1=strain[3], u2=UNSET, u3=UNSET, amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_22', createStepName='Step-1',
                                region=assembly.sets['Ref-Y'], u1=UNSET, u2=strain[1], u3=UNSET, amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_23', createStepName='Step-1',
                                region=assembly.sets['Ref-Y'], u1=UNSET, u2=UNSET, u3=strain[4], amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        
        #### strain for the third dummy point(the top one E31 E32 E33)
        model.DisplacementBC(name='E_31', createStepName='Step-1',
                                region=assembly.sets['Ref-Z'], u1=strain[5], u2=UNSET, u3=UNSET, amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_32', createStepName='Step-1',
                                region=assembly.sets['Ref-Z'], u1=UNSET, u2=strain[4], u3=UNSET, amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_33', createStepName='Step-1',
                                region=assembly.sets['Ref-Z'], u1=UNSET, u2=UNSET, u3=strain[2], amplitude=UNSET, fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
    else:
        # create strain path 
        num_point = strain_amplitude.shape[0]
        names = ['E_11', 'E_22', 'E_33', 'E_12', 'E_23', 'E_13']
        for ii, name in enumerate(names):
            path_table_temp = np.zeros((num_point, 2))
            path_table_temp[:, 0] = np.linspace(
                0, time_period, num_point, endpoint=True
            )
            path_table_temp[:, 1] = strain_amplitude[:, ii]
            create_random_path(model=model, path_table=path_table_temp, path_name=name)
        
        # apply the strain path to the model
        #### strain for the first dummy point(the front one E11,E12,E13)
        model.DisplacementBC(name='E_11', createStepName='Step-1',
                                region=assembly.sets['Ref-X'], u1=strain[0], u2=UNSET, u3=UNSET, amplitude='E_11', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_12', createStepName='Step-1',
                                region=assembly.sets['Ref-X'], u1=UNSET, u2=strain[3], u3=UNSET, amplitude='E_12', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_13', createStepName='Step-1',
                                region=assembly.sets['Ref-X'], u1=UNSET, u2=UNSET, u3=strain[5], amplitude='E_13', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        
        #### strain for the second dummy point(the front one E21,E22,E23)
        model.DisplacementBC(name='E_21', createStepName='Step-1',
                                region=assembly.sets['Ref-Y'], u1=strain[3], u2=UNSET, u3=UNSET, amplitude='E_12', fixed=OFF,   
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_22', createStepName='Step-1',
                                region=assembly.sets['Ref-Y'], u1=UNSET, u2=strain[1], u3=UNSET, amplitude='E_22', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_23', createStepName='Step-1',
                                region=assembly.sets['Ref-Y'], u1=UNSET, u2=UNSET, u3=strain[4], amplitude='E_23', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)

        #### strain for the third dummy point(the top one E31 E32 E33)
        model.DisplacementBC(name='E_31', createStepName='Step-1',
                                region=assembly.sets['Ref-Z'], u1=strain[5], u2=UNSET, u3=UNSET, amplitude='E_13', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_32', createStepName='Step-1',
                                region=assembly.sets['Ref-Z'], u1=UNSET, u2=strain[4], u3=UNSET, amplitude='E_23', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
        model.DisplacementBC(name='E_33', createStepName='Step-1',
                                region=assembly.sets['Ref-Z'], u1=UNSET, u2=UNSET, u3=strain[2], amplitude='E_33', fixed=OFF,
                                distributionType=UNIFORM, fieldName='', localCsys=None)
    # restrict the rigid movement
    assembly.SetFromNodeLabels(
        name="supportnode",
        nodeLabels=(("RVE", (support[0].label,)),),
        unsorted=True,
    )

    model.DisplacementBC(
        amplitude=UNSET,
        createStepName="Step-1",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="rigid",
        region=assembly.sets["supportnode"],
        u1=0.0,
        u2=0.0,
        u3=0.0,
    )
    # create job
    mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=num_cpus, numDomains=num_cpus, numGPUs=0)
    mdb.jobs[job_name].writeInput(consistencyChecking=OFF)

def post_process(dict):
    job_name = dict['job_name']
    
    # size of the rve 
    size = dict['len_end']- dict['len_start'] - 2*dict['radius_mu']

    rve_odb = openOdb(path=job_name+ '.odb')

    # get the element set
    entire_element_set = rve_odb.rootAssembly.elementSets[" ALL ELEMENTS"]
    ref_x_node_set = rve_odb.rootAssembly.nodeSets["REF-X"]
    ref_y_node_set = rve_odb.rootAssembly.nodeSets["REF-Y"]
    ref_z_node_set = rve_odb.rootAssembly.nodeSets["REF-Z"]

    # volume of the RVE 
    volume = size**3
    # get the steps
    my_steps = rve_odb.steps

    total_frames = 0
    for ii in range(len(my_steps)):
        total_frames = total_frames +  len(my_steps[my_steps.keys()[ii]].frames)

    rve_frame = rve_odb.steps[my_steps.keys()[0]].frames[0]


    # define variables
    strain = np.zeros((total_frames, 3, 3))
    stress = np.zeros((total_frames, 3, 3))
    total_time = np.zeros(total_frames)

    # displacement
    U_ref_x = np.zeros((total_frames, 3))
    U_ref_y = np.zeros((total_frames, 3))
    U_ref_z = np.zeros((total_frames, 3))

    # reaction forces
    RF_ref_x = np.zeros((total_frames, 3))
    RF_ref_y = np.zeros((total_frames, 3))
    RF_ref_z = np.zeros((total_frames, 3))

    # plastic energy 

    # loop over all frames
    for ii in range(total_frames):
        rve_frame = rve_odb.steps[my_steps.keys()[0]].frames[ii]
        total_time[ii] = rve_frame.frameValue
        # get the field output
        fo = rve_frame.fieldOutputs
        # get the displacement
        U_ref_x[ii] = np.array(fo['U'].getSubset(region=ref_x_node_set).values[0].data)
        U_ref_y[ii] = np.array(fo['U'].getSubset(region=ref_y_node_set).values[0].data)
        U_ref_z[ii] = np.array(fo['U'].getSubset(region=ref_z_node_set).values[0].data)
        # get the reaction forces
        RF_ref_x[ii] = np.array(fo['RF'].getSubset(region=ref_x_node_set).values[0].data)
        RF_ref_y[ii] = np.array(fo['RF'].getSubset(region=ref_y_node_set).values[0].data)
        RF_ref_z[ii] = np.array(fo['RF'].getSubset(region=ref_z_node_set).values[0].data)

        # post process the field output to get the homogenized strain and stress 
        # xx, xy, xz
        # yx, yy, yz
        # zx, zy, zz
        strain[ii,0,: ] = U_ref_x[ii]
        strain[ii,1,: ] = U_ref_y[ii]
        strain[ii,2,: ] = U_ref_z[ii]

        # stress
        stress[ii,0,:] = RF_ref_x[ii]/volume
        stress[ii,1,:] = RF_ref_y[ii]/volume
        stress[ii,2,:] = RF_ref_z[ii]/volume


    results = { 'strain': strain,
                'stress': stress,
                'total_time': total_time}

    with open('results.pkl', 'wb') as f:
        pickle.dump(results, f, )