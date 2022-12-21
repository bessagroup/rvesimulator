# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *

# import other class


class AbaqusMaterialLib:

    def __init__(self, name_mat, model, part, name_set):
        self.model = model
        self.part = part
        self.name_mat = name_mat
        self.name_set = name_set

    def CreateElasticMaterial(self, E, v):
        material = self.model.Material(name=self.name_mat)
        material.Elastic(table=((E, v),))
        self.model.HomogeneousSolidSection(name=self.name_mat, material=self.name_mat, thickness=None)
        self.part.SectionAssignment(region=self.part.sets[self.name_set], sectionName=self.name_mat, offset=0.0,
                                    offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

        return material

    def CreateVonMisesPlasticMaterial(self, E, v, yield_criterion):
        material = self.model.Material(name=self.name_mat)
        material.Elastic(table=((E, v),))
        material.Plastic(table=(yield_criterion))
        self.model.HomogeneousSolidSection(name=self.name_mat, material=self.name_mat, thickness=None)
        self.part.SectionAssignment(region=self.part.sets[self.name_set], sectionName=self.name_mat, offset=0.0,
                                    offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def CreateUserMaterialMirkhalafModel(self, density, num_stv, para):
        """
        :param density: density of the polymer
        :param num_stv: number of state variables
        :param para:  list of parameters for Mirkhalaf model, 12 paramters in total
        :return: nothing but will update the model.material and assign material property to
        corresponding part 
        """

        material = self.model.Material(name=self.name_mat)
        material.Density(table=((density,),))
        material.Depvar(n=num_stv)
        material.UserMaterial(mechanicalConstants=(para[0], para[1], para[2], para[3], para[4], para[5],
                                                   para[6], para[7], para[8], para[9], para[10], para[11]))
        self.model.HomogeneousSolidSection(name=self.name_mat, material=self.name_mat, thickness=None)
        self.part.SectionAssignment(region=self.part.sets[self.name_set], sectionName=self.name_mat, offset=0.0,
                                    offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

