# abaqus
from json import loads
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
        
        if E22 == E12 == 0.0 : 
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



