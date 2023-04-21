# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *

# import other class


class MaterialLib:
    def create_elastic_material(
        self, geometry_set_name, material_name, youngs, poisson
    ):
        """create elastic material 

        Parameters
        ----------
        geometry_set_name : str
            geometry set name 
        material_name : str
            material name 
        youngs : float
            youngs modulus
        poisson : float
            poisson ratio
        """

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


    def create_user_material_VEVP_Leonov_model(self,
                                                geometry_set_name,
                                                material_name,
                                                num_stv,
                                                para):
        """use VEVP Leonov model as user material 

        Parameters
        ----------
        geometry_set_name : _type_
            _description_
        material_name : _type_
            _description_
        num_stv : _type_
            _description_
        para : _type_
            _description_
        """
        material = self.model.Material(name=material_name)
        material.Depvar(n=num_stv)
        material.UserMaterial(
            mechanicalConstants=(
                para[0],
                para[1],
                para[2],
                para[3],
                para[4],
                para[5],
                para[6],
                para[7],
                para[8],
                para[9],
                para[10],
                para[11],
                para[12],
                para[13],
            )
        )
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


    def create_user_material_VP_Leonov_model(self,
                                            geometry_set_name,
                                            material_name,  
                                            num_stv, 
                                            para):
        """use VP Leonov model as user material 

        Parameters
        ----------
        geometry_set_name : _type_
            _description_
        material_name : _type_
            _description_
        num_stv : _type_
            _description_
        para : _type_
            _description_
        """
        material = self.model.Material(name=material_name)
        material.Depvar(n=num_stv)
        material.UserMaterial(
            mechanicalConstants=(
                para[0],
                para[1],
                para[2],
                para[3],
                para[4],
                para[5],
                para[6],
                para[7],
                para[8],
                para[9],
                para[10],
            )
        )
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