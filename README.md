<p align="center">
    <img src="docs/images/logo.png" width="70%" align="center">
</p>

## rvesimulator

RVE simulation via ABAQUS[ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/) secondary development.

---

[![Python](https://img.shields.io/pypi/pyversions/f3dasm)](https://www.python.org)
[![GitHub license](https://img.shields.io/badge/license-BSD-blue)](https://github.com/bessagroup/rvesimulator)

[**Docs**](https://bessagroup.github.io/rvesimulator/)
| [**Installation**](https://bessagroup.github.io/rvesimulator/get_started.html)
| [**GitHub**](https://github.com/bessagroup/rvesimulator)

---

### **Summary**

This repo aims to do representative volume element (RVE) simulation via [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/). It has two main parts, the first part is to generate microstructure of representative volume element (RVE), and the second part is realize the functionality of do RVE simulation automatically and has the ability to generate large dataset without any GUI operation.

This idea of the repo is shown by the following figures:

<p align="center">
    <img src="docs/images/rvesimulator.svg" width="70%" align="center">
</p>

---

### **State of need**

In order to use this repo, one needs first to know basics of ABAQUS secondary development, and basic under standing of design of experiment (DoE).

---

### **Authorship**

**Author**:

- Jiaxiang Yi ([J.Yi@tudelft.nl](mailto:J.Yi@tudelft.nl))

**Author filliation**:

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

(3). install dependancies

```
pip install -r requirements.txt
```

**Contents**

Part I: microstructure generation

Part II: rve simulation

---

### **Community Support**

Now everything is under developing, if you have any question, please contact to the developer.

### **Credits**

### **License**

Copyright 2020, Jiaxiang Yi

All rights reserved.

rvesimulator is a free and open-source software published under a MIT License.
