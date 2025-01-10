"""rvesimulator - python package for simulating Representative Volume Elements
(RVEs) via commercial finite element software Abaqus.

This package provides tools for generating micro-structures, simulating RVEs,
and post-processing simulation results.

Usage
-----

>>> import rvesimulator

Links
-----

- Documentation: https://bessagroup.github.io/rvesimulator/index.html

Author: Jiaxiang Yi (J.Yi@tudelft.nl)
"""

#                                                        Authorship and Credits
# =============================================================================
__author__ = 'Jiaxiang Yi (J.Yi@tudelft.nl)'
__credits__ = ['Jiaxiang Yi']
__status__ = 'Stable'
#
# =============================================================================
from .abaqus2py import AbaqusSimulator
from .additions import ampitudesampler, hardening_law
from .benchmarks.asca_rve import ASCA_RVE
from .benchmarks.cddm_rve import CDDM_RVE
from .benchmarks.hollow_plate_sve import (ElasticRegularLoads,
                                          VonMisesPlasticPathLoads,
                                          VonMisesPlasticRegularLoads)
from .microstructure import circle_particles, sphere_particles
