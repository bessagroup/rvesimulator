# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *


class HollowPlate:
    def create_part(self):
        # create sketch 
        sketch = self.model.ConstrainedSketch(
            name="sketch_profile", sheetSize=2*self.length
        )
        # create rectangle  
        sketch.rectangle(point1=(-self.length / 2, -self.width / 2), 
                         point2=(self.length / 2, self.width / 2))
        # create the circle inside 
        sketch.CircleByCenterPerimeter(
            center=self.center, point1=(self.radius, 0.0)
        )
        # create part 
        self.part = self.model.Part(
            name=self.part_name,
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part.BaseShell(sketch=sketch)


class MultiCirclesPlates: 
    def create_part(self): 
        # create sketch 
        sketch = self.model.ConstrainedSketch(
            name="sketch_profile", sheetSize=2.0
        )
        # create the rectangle 
        sketch.rectangle(
            point1=(0, 0),
            point2=(self.length, self.width),
        )
        # create circles inside the rectangle 
        for ii in range(len(self.circles_information)):
            sketch.CircleByCenterPerimeter(
                center=(self.circles_information[ii][0], self.circles_information[ii][1]),
                point1=(
                    self.circles_information[ii][0] + self.circles_information[ii][2],
                    self.circles_information[ii][1],
                ),
            )
        # create part based on the geometry
        self.part = self.model.Part(
            name=self.part_name,
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part.BaseShell(sketch=sketch)         


# =============================================================================
class MultiCirclesFullyInclusion:

    """statistical rve, all circle fibers are inside the square, 2 materials"""

    def create_part(self):

        # Create the  fibers
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        for ii in range(len(self.circles_information)):
            sketch.CircleByCenterPerimeter(
                center=(self.circles_information[ii][0], self.circles_information[ii][1]),
                point1=(self.circles_information[ii][0] + self.circles_information[ii][2], 
                        self.circles_information[ii][1]))
            
        self.model.ConstrainedSketch(
            name="Fibre_Sketch_Prov", objectToCopy=sketch
        )
        self.model.sketches.changeKey(
            fromName="__profile__", toName="Fibre_Sketch_Prov"
        )
        # Matrix Sketch
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.rectangle(
            point1=(0.0, 0.0),
            point2=(self.length, self.width),
        )
        self.model.sketches.changeKey(
            fromName="__profile__", toName="Matrix_Sketch"
        )

        # fibers Part
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.retrieveSketch(sketch=self.model.sketches["Fibre_Sketch_Prov"])

        self.part = self.model.Part(
            name="Fibre_Final_Part",
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part = self.model.parts["Fibre_Final_Part"]
        self.part.BaseShell(sketch=sketch)
        self.part = self.model.parts["Fibre_Final_Part"]

        del self.model.sketches["__profile__"]
        del self.model.sketches["Fibre_Sketch_Prov"]

        # The original square part (matrix part)
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.retrieveSketch(sketch=self.model.sketches["Matrix_Sketch"])
        self.part = self.model.Part(
            name="Matrix_Part",
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part = self.model.parts["Matrix_Part"]
        self.part.BaseShell(sketch=sketch)
        self.part = self.model.parts["Matrix_Part"]
        del self.model.sketches["__profile__"]
        del self.model.sketches["Matrix_Sketch"]

        # Fibre Part
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Fibre_Final_Part"]
        self.assembly.Instance(
            name="Fibre_Final_Instance", part=self.part, dependent=ON
        )
        # Matrix Part
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Matrix_Part"]
        self.assembly.Instance(
            name="Matrix_Instance", part=self.part, dependent=ON
        )
        self.assembly.InstanceFromBooleanCut(
            name="Matrix_Final_Part",
            instanceToBeCut=self.assembly.instances["Matrix_Instance"],
            cuttingInstances=(
                self.assembly.instances["Fibre_Final_Instance"],
            ),
            originalInstances=DELETE,
        )
        self.assembly.makeIndependent(
            instances=(self.assembly.instances["Matrix_Final_Part-1"],)
        )
        self.assembly.features.changeKey(
            fromName="Matrix_Final_Part-1", toName="Matrix_Final_Instance"
        )
        del self.model.parts["Matrix_Part"]

        # Assembly
        self.part = self.model.parts["Fibre_Final_Part"]
        self.assembly.Instance(
            name="Fibre_Final_Instance", part=self.part, dependent=ON
        )
        self.assembly.InstanceFromBooleanMerge(
            name="Final_Stuff",
            instances=(
                self.assembly.instances["Matrix_Final_Instance"],
                self.assembly.instances["Fibre_Final_Instance"],
            ),
            keepIntersections=ON,
            originalInstances=DELETE,
            domain=GEOMETRY,
        )
        del self.model.parts["Matrix_Final_Part"]
        del self.model.parts["Fibre_Final_Part"]


class MultiCirclesInclusion: 

    def create_part(self):
        """create part

        Returns
        -------
        abaqus part
            abaqus part
        """

        # Create the  fibers
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        for ii in range(len(self.circles_information)):
            sketch.CircleByCenterPerimeter(
                center=(self.circles_information[ii][0], self.circles_information[ii][1]),
                point1=(self.circles_information[ii][0] + self.circles_information[ii][2], self.circles_information[ii][1]),
            )
        self.model.ConstrainedSketch(
            name="Fibre_Sketch_Prov", objectToCopy=sketch
        )
        self.model.sketches.changeKey(
            fromName="__profile__", toName="Fibre_Sketch_Prov"
        )

        # create a frame to cut the fibers so they can stay inside the square
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.rectangle(
            point1=(self.len_start + self.radius_mu, self.wid_start + self.radius_mu),
            point2=(self.len_end - self.radius_mu, self.wid_end - self.radius_mu),
        )
        sketch.rectangle(
            point1=(
                self.len_start - 2 * self.radius_mu,
                self.wid_start - 2 * self.radius_mu,
            ),
            point2=(
                self.len_end + 2 * self.radius_mu,
                self.wid_end + 2 * self.radius_mu,
            ),
        )
        self.model.sketches.changeKey(
            fromName="__profile__", toName="Fibre_Sketch_Trim"
        )

        #  Matrix Sketch
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.rectangle(
            point1=(self.len_start + self.radius_mu, self.wid_start + self.radius_mu),
            point2=(self.len_end - self.radius_mu, self.wid_end - self.radius_mu),
        )
        self.model.sketches.changeKey(
            fromName="__profile__", toName="Matrix_Sketch"
        )

        # Fibres Part (Provisional): before trim
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.retrieveSketch(sketch=self.model.sketches["Fibre_Sketch_Prov"])
        self.part = self.model.Part(
            name="Fibre_Part_Prov",
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part = self.model.parts["Fibre_Part_Prov"]
        self.part.BaseShell(sketch=sketch)
        self.part = self.model.parts["Fibre_Part_Prov"]
        del self.model.sketches["__profile__"]
        del self.model.sketches["Fibre_Sketch_Prov"]

        # trim part
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.retrieveSketch(sketch=self.model.sketches["Fibre_Sketch_Trim"])
        self.part = self.model.Part(
            name="Fibre_Part_Trim",
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part = self.model.parts["Fibre_Part_Trim"]
        self.part.BaseShell(sketch=sketch)
        self.part = self.model.parts["Fibre_Part_Trim"]
        del self.model.sketches["__profile__"]
        del self.model.sketches["Fibre_Sketch_Trim"]

        # The original square part (matrix part)
        sketch = self.model.ConstrainedSketch(
            name="__profile__", sheetSize=self.length
        )
        sketch.retrieveSketch(sketch=self.model.sketches["Matrix_Sketch"])
        self.part = self.model.Part(
            name="Matrix_Part",
            dimensionality=TWO_D_PLANAR,
            type=DEFORMABLE_BODY,
        )
        self.part = self.model.parts["Matrix_Part"]
        self.part.BaseShell(sketch=sketch)
        self.part = self.model.parts["Matrix_Part"]
        del self.model.sketches["__profile__"]
        del self.model.sketches["Matrix_Sketch"]

        # Fibre Part
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Fibre_Part_Prov"]
        self.assembly.Instance(
            name="Fibre_Instance_Prov", part=self.part, dependent=ON
        )
        self.part = self.model.parts["Fibre_Part_Trim"]
        self.assembly.Instance(
            name="Fibre_Instance_Trim", part=self.part, dependent=ON
        )
        self.assembly.InstanceFromBooleanCut(
            name="Fibre_Final_Part",
            instanceToBeCut=self.assembly.instances["Fibre_Instance_Prov"],
            cuttingInstances=(self.assembly.instances["Fibre_Instance_Trim"],),
            originalInstances=DELETE,
        )
        self.assembly.makeIndependent(
            instances=(self.assembly.instances["Fibre_Final_Part-1"],)
        )
        self.model.rootAssembly.features.changeKey(
            fromName="Fibre_Final_Part-1", toName="Fibre_Final_Instance"
        )
        del self.model.parts["Fibre_Part_Prov"]
        del self.model.parts["Fibre_Part_Trim"]

        # Matrix Part
        self.assembly = self.model.rootAssembly
        self.part = self.model.parts["Matrix_Part"]
        self.assembly.Instance(
            name="Matrix_Instance", part=self.part, dependent=ON
        )
        self.assembly.InstanceFromBooleanCut(
            name="Matrix_Final_Part",
            instanceToBeCut=self.assembly.instances["Matrix_Instance"],
            cuttingInstances=(
                self.assembly.instances["Fibre_Final_Instance"],
            ),
            originalInstances=DELETE,
        )
        self.assembly.makeIndependent(
            instances=(self.assembly.instances["Matrix_Final_Part-1"],)
        )
        self.assembly.features.changeKey(
            fromName="Matrix_Final_Part-1", toName="Matrix_Final_Instance"
        )
        del self.model.parts["Matrix_Part"]

        # Assembly and Meshing
        self.part = self.model.parts["Fibre_Final_Part"]
        self.assembly.Instance(
            name="Fibre_Final_Instance", part=self.part, dependent=ON
        )
        self.assembly.InstanceFromBooleanMerge(
            name="Final_Stuff",
            instances=(
                self.assembly.instances["Matrix_Final_Instance"],
                self.assembly.instances["Fibre_Final_Instance"],
            ),
            keepIntersections=ON,
            originalInstances=DELETE,
            domain=GEOMETRY,
        )
        del self.model.parts["Matrix_Final_Part"]
        del self.model.parts["Fibre_Final_Part"]
