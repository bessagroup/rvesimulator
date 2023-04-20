# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *

# import other class


class MaterialLib:
    def create_elastic_material(
        self, geometry_set_name, material_name, youngs, poisson
    ):

        self.material = self.model.Material(name=material_name)
        self.material.Elastic(table=((youngs, poisson),))
        self.model.HomogeneousSolidSection(
            name=material_name, material=material_name, thickness=None
        )
        self.part.SectionAssignment(
            region=self.part.sets[geometry_set_name],
            sectionName=material_name,
            offset=0.0,
            offsetType=MIDDLE_SURFACE,
            offsetField="",
            thicknessAssignment=FROM_SECTION,
        )

    def create_von_mises_plastic_material(self,
                                          geometry_set_name,
                                          material_name,
                                          youngs,
                                          poisson,
                                          hardening_table):

        material = self.model.Material(name=material_name)
        material.Elastic(table=((youngs, poisson),))
        material.Plastic(table=(hardening_table))
        self.model.HomogeneousSolidSection(
            name=material_name, material=material_name, thickness=None
        )
        self.part.SectionAssignment(
            region=self.part.sets[geometry_set_name],
            sectionName=material_name,
            offset=0.0,
            offsetType=MIDDLE_SURFACE,
            offsetField="",
            thicknessAssignment=FROM_SECTION,
        )
