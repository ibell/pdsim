# PDSim

What is PDSim?

* A python/cython-based library for Positive-Displacement machine Simulation
* Written in a flexible way so that a wide variety of machines can be implemented with the same core. 
* Two papers about PDSim are available in *Int. J. Refrig.*:
    * [PDSim: A general quasi-steady modeling approach for positive displacement compressors and expanders](https://doi.org/10.1016/j.ijrefrig.2019.09.002) preprint: [preprints/PaperI.pdf](doc/papers/PaperI.pdf)

This version uses CoolProp v6.x

## Documentation/Examples

* Try me on binder [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ibell/pdsim/master?filepath=doc%2Fnotebooks)

* Build status on TravisCI: [![Build Status](https://travis-ci.com/ibell/pdsim.svg?branch=master)](https://travis-ci.com/ibell/pdsim)

* The documention is available on ReadTheDocs: [![Documentation Status](https://readthedocs.org/projects/pdsim/badge/?version=latest)](http://pdsim.readthedocs.io/en/latest/?badge=latest)

## Install

Requirements:
* numpy
* cython
* matplotlib

To setup a conda environment, you can do move into the :
```
conda env create
```

```
pip -vvv install .
```
