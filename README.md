<p align="center">
    <img src="docs/images/logo.png" width="70%" align="center">
</p>

## rvesimulator

RVE simulation via [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/) secondary development.

---

[**Documentation**](https://bessagroup.github.io/rvesimulator/)
| [**Installation**](https://bessagroup.github.io/rvesimulator/get_started.html)
| [**GitHub**](https://github.com/bessagroup/rvesimulator)
| [**Tutorials**](https://github.com/bessagroup/rvesimulator/tree/main/tutorials)

---

### **Summary**

The **rvesimulator** aims to provide a user-friendly, automated Python-based framework conducting Representative Volume Element (RVE) simulation via powerful Finite Element Method (FEM) software Abaqus. By utilizing this repository, large amount of reliable FEM data-set generation is possible with RVEs encompassing materials from elastic to plastic composites .

**rvesimulator** provides:

1.  A cross-platform function to run arbitrary Python-Abaqus script without graphical user interface (GUI), it offers users a convenience way to run their unique scripts
2.  Python-Abaqus scripts to simulate RVE with different design of experiments including various micro-structures, material laws, and loading;
3.  Benchmarks of running prevalent RVEs covering elastic, hyper-elastic, plastic materials are provided, which illustrates the general pipeline (preprocess, execution, and postprocess) of the developed framework

This idea of the repo is shown by the following figures:

<p align="center">
    <img src="docs/images/rvesimulator.svg" width="70%" align="center">
</p>

---

### **State of need**

In order to use this repo, one needs first to know basics of ABAQUS secondary development, and basic under standing of design of experiment (DoE) so that one can design their own RVEs.

---

### **Authorship**

**Author**:

- Jiaxiang Yi ([J.Yi@tudelft.nl](mailto:J.Yi@tudelft.nl))

**Author affiliation**:

- Delft University of Technology

---

### Get started

**Installation**

(1). git clone the repo to your local machine

```
git clone https://github.com/bessagroup/rvesimulator.git
```

(2). go to the local folder where you cloned the repo, and pip install it with editable mode

```
pip install --verbose --no-build-isolation --editable .
```

(3). install dependencies

```
pip install -r requirements.txt
```

<!-- **Contents**

Part I: microstructure generation

Part II: rve simulation -->

---

### **Community Support**

Now everything is under developing, if you have any question, please contact to the developer.

### **Credits**

### **License**

Copyright 2020, Jiaxiang Yi

All rights reserved.

rvesimulator is a free and open-source software published under a MIT License.
