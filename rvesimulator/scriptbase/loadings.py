# abaqus
import json

import numpy
from abaqus import *
from abaqusConstants import *
from caeModules import *


class Loading2D:

    def __init__(self, model, assembly, instance_name):
        self.model = model
        self.assembly = assembly
        self.instance = instance_name 
         
    def create_loads(self, loads):
        # loads should be a list with three elemnts 
        E11 = loads[0] 
        E22 = loads[1]
        E12 = loads[2]      
        if E22 == E12 == 0.0: 
            # x axis tension 
            self.E11(E11=E11)
            self.RigidConsE11() 
        elif E11 == E12 == 0.0: 
            # y axis tension 
            self.E22(E22=E22) 
            self.RigidConsE22()
        elif E11 == E22 ==0.0: 
            # pure shear 
            self.E12(E12=E12) 
        else: 
            self.Complex(loads=loads) 
            
    def E11(self, E11):
        self.model.DisplacementBC(name='E_11', createStepName='Step-1',
                                  region=self.assembly.sets['Ref-R'], u1=E11, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,
                                  distributionType=UNIFORM, fieldName='', localCsys=None)
    def E22(self, E22):
        self.model.DisplacementBC(name='E_22', createStepName='Step-1',
                                  region=self.assembly.sets['Ref-T'], u1=UNSET, u2=E22, ur3=UNSET, amplitude=UNSET, fixed=OFF,
                                  distributionType=UNIFORM, fieldName='', localCsys=None)
    def E12(self, E12):
        self.model.DisplacementBC(name='E_12', createStepName='Step-1',
                                 region=self.assembly.sets['Ref-R'], u1=UNSET, u2=E12, ur3=UNSET, amplitude=UNSET, fixed=OFF,
                                 distributionType=UNIFORM, fieldName='', localCsys=None)
        self.model.DisplacementBC(name='E_21', createStepName='Step-1',
                                 region=self.assembly.sets['Ref-T'], u1=E12, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,
                                 distributionType=UNIFORM, fieldName='', localCsys=None)
    def Complex(self, loads): 
        self.model.DisplacementBC(name='E_R_Ref', createStepName='Step-1',
                                  region=self.assembly.sets['Ref-R'], u1=loads[0], u2=loads[2], ur3=UNSET, amplitude=UNSET, fixed=OFF,
                                  distributionType=UNIFORM, fieldName='', localCsys=None) 
        self.model.DisplacementBC(name='E_T_Ref', createStepName='Step-1',
                            region=self.assembly.sets['Ref-T'], u1=loads[2], u2=loads[1], ur3=UNSET, amplitude=UNSET, fixed=OFF,
                            distributionType=UNIFORM, fieldName='', localCsys=None) 
            
    def RigidConsE11(self):
        self.model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='',
                                  fixed=OFF, localCsys=None, name='RigidConsE111', region=self.assembly.instances[self.instance].sets['VertexLB'], u1=0.0, u2=UNSET, ur3=UNSET)
        self.model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='',
                                  fixed=OFF, localCsys=None, name='RigidConsE112', region=self.assembly.instances[self.instance].sets['VertexLT'], u1=0.0, u2=0.0, ur3=UNSET)

    def RigidConsE22(self):
        self.model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='',
                                  fixed=OFF, localCsys=None, name='RigidConsE221', region=self.assembly.instances[self.instance].sets['VertexLB'], u1=0.0, u2=0.0, ur3=UNSET)
        self.model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='',
                                  fixed=OFF, localCsys=None, name='RigidConsE222', region=self.assembly.instances[self.instance].sets['VertexRB'], u1=UNSET, u2=0.0, ur3=UNSET)
    
    def create_path_load(self, load_magnitude, path_table, time_period):
        # loads should be a list with three elemnts
        E11 = load_magnitude[0]
        E22 = load_magnitude[1]
        E12 = load_magnitude[2]
        num_point = path_table.shape[0]

        # create the path amplitude first 
        path_name = ['xx','yy','xy'] 
        for ii, name in enumerate(path_name):
            path_table_temp = numpy.zeros((num_point, 2))
            path_table_temp[:,0] = numpy.linspace(0, time_period, num_point, endpoint=True)
            path_table_temp[:,1] = path_table[:, ii]
            self._create_random_path(path_table_temp, name)  
    
        # create load for corresonding part 
        # create E_11
        self.model.DisplacementBC(name='E11', createStepName='Step-1', 
                    region=self.assembly.sets["Ref-R"], u1=E11, u2=UNSET, ur3=UNSET, amplitude=path_name[0], fixed=OFF, 
                    distributionType=UNIFORM, fieldName='', localCsys=None) 
        # create E_22 
        self.model.DisplacementBC(name='E22', createStepName='Step-1', 
                    region=self.assembly.sets["Ref-T"], u1=UNSET, u2=E22, ur3=UNSET, amplitude=path_name[1], fixed=OFF, 
                    distributionType=UNIFORM, fieldName='', localCsys=None) 
        # create E_12
        self.model.DisplacementBC(name='E12', createStepName='Step-1', 
                    region=self.assembly.sets["Ref-R"], u1=UNSET, u2=E12, ur3=UNSET, amplitude=path_name[2], fixed=OFF, 
                    distributionType=UNIFORM, fieldName='', localCsys=None) 
        # create E_21
        self.model.DisplacementBC(name='E21', createStepName='Step-1', 
                    region=self.assembly.sets["Ref-T"], u1=E12, u2=UNSET, ur3=UNSET, amplitude=path_name[2], fixed=OFF, 
                    distributionType=UNIFORM, fieldName='', localCsys=None)     

    def _create_random_path(self, path_tale, path_name):

        # generate the table
        self.model.TabularAmplitude(
            name=path_name,
            timeSpan=TOTAL,
            smooth=SOLVER_DEFAULT,
            data=(path_tale))