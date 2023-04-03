## abaqus rve simulator 

This repo aims to do representative volume element (RVE) simulation via [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/). It has two main parts, the first part is to generate microstructure of representative volume element (RVE), and the second part is realize the functionality of do RVE simulation automatically and has the ability to generate large dataset without any GUI operation. 

---
**Author**:
- Jiaxiang Yi([J.Yi@tudelft.nl](mailto:J.Yi@tudelft.nl)) 

**Author filliation**:
- Delft University of Technology 
---
### 1. Installation 
1. git clone the repo to your local machine 
``` 
git clone https://github.com/bessagroup/rvesimulator.git 
```
2. go to the local folder where you cloned the repo, and pip install it with editable mode 
```
pip install --verbose --no-build-isolation --editable .
```
3. install dependancies 
```
pip install -r requirements.txt
```
**Note**: this repo also dependent on [f3dasm](https://github.com/bessagroup/F3DASM/tree/main), but you have to install f3dasm manually according to its guidance 

---
### 2. Contents
#### Part I: microstructure generation

```from abaqusrve.microstructure.circle_particles import CircleParticles
microstructure_generator = CircleParticles(
    length=1.0,
    width=1.0,
    radius_mu=0.1,
    radius_std=0.03,
    vol_req=0.45,
)
microstructure_generator.generate_microstructure()
```
 
#### Part I: rve simulation
for more details of this part, please go to the tutorials folder to see details instructions for different well defined rve simulation cases. 