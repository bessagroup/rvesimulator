<p align="center">
    <img src="docs/images/logo.png" width="70%" align="center">
</p>

---

[**Documentation**](https://bessagroup.github.io/rvesimulator/)
| [**Installation**](https://bessagroup.github.io/rvesimulator/get_started.html)
| [**GitHub**](https://github.com/bessagroup/rvesimulator)
| [**Tutorials**](https://github.com/bessagroup/rvesimulator/tree/main/tutorials)

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

[rvesimulator](https://github.com/bessagroup/rvesimulator) is an open-source Python-based framework for RVE simulation via [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/) secondary development. It is noted that ABAQUS 2021 is used in the development. In order to conduct RVE simulation via [rvesimulator](https://github.com/bessagroup/rvesimulator) , basic knowledge of [ABAQUS](https://www.3ds.com/products-services/simulia/products/abaqus/) usage, especially Python scripting, is required. Meanwhile, basic understanding of RVE simulation is also required. Several examples are provided in the [Tutorials](https://github.com/bessagroup/rvesimulator/tree/main/tutorials) where users can learn the basic usage of [rvesimulator](https://github.com/bessagroup/rvesimulator). If the user wants to develop their own RVE simulation, the [Documentation](https://bessagroup.github.io/rvesimulator/) provides a detailed description on how to adapt the framework to their own needs.

Consider leaving a star if you think [rvesimulator](https://github.com/bessagroup/rvesimulator) is useful for the research community!

---

### **Authorship & Citation**

**Author**:

- Jiaxiang Yi ([J.Yi@tudelft.nl](mailto:J.Yi@tudelft.nl))

**Author affiliation**:

- Delft University of Technology

---

If you use [rvesimulator](https://github.com/bessagroup/rvesimulator), please cite the following [paper](https://openreview.net/forum?id=511z1DGjPi):

```
    @inproceedings{
    yi2023rvesimulator,
    title={rvesimulator: An automated representative volume element simulator for data-driven material discovery},
    author={Jiaxiang Yi and Miguel Anibal Bessa},
    booktitle={AI for Accelerated Materials Design - NeurIPS 2023 Workshop},
    year={2023},
    url={https://openreview.net/forum?id=511z1DGjPi}
    }
```

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

Now the [rvesimulator](https://github.com/bessagroup/rvesimulator) is still under development, if you have any questions, please feel free to contact [Jiaxiang Yi](mailto:J.Yi@tudelf.nl), or open an issue on the [GitHub](https://github.com/bessagroup/rvesimulator/issues) issues page.

<!-- ### **Credits** -->

### **License**

Copyright 2022, Jiaxiang Yi

All rights reserved.

rvesimulator is a free and open-source software published under a MIT License.
