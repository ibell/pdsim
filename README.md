# PDSim

What is PDSim?

* A python/cython-based library for *P*ositive-*D*isplacement machine *Sim*ulation
* Written in a flexible way so that a wide variety of machines can be implemented with the same core. 
* Two papers about PDSim are available in *Int. J. Refrig.* in 2020:
    * [PDSim: A general quasi-steady modeling approach for positive displacement compressors and expanders](https://doi.org/10.1016/j.ijrefrig.2019.09.002) preprint: [PaperI](doc/papers/PaperI.pdf)
    * [PDSim: A general quasi-steady modeling approach for positive displacement compressors and expanders](https://doi.org/10.1016/j.ijrefrig.2019.10.015) preprint: [PaperII](doc/papers/PaperII.pdf)

This version uses CoolProp v6.x

# Credits/History

The core code of PDSim was written by Ian Bell, as part of a collaborative project with the [Ray W. Herrick Lab faculty](https://engineering.purdue.edu/Herrick) while backpacking around the world. It has its origins in the work done on [modeling scroll compressors](https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1003), and to this day the most comprehensive analysis is possible for scroll compressors.  Examples are available of:

* [hermetic scroll compressor](examples/scroll_compressor.py)
* [hermetic scroll compressor w/ vapor injection](examples/scroll_compressor_w_VI.py)
* [hermetic scroll compressor w/ discharge valve](examples/scroll_compressor_valve.py)

Examples are also available of:

* [reciprocating compressor](examples/recip_compressor.py)
* [piston expander](examples/piston_expander.py)

## Documentation/Examples

* Try me on binder [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ibell/pdsim/master?filepath=doc%2Fnotebooks)

* Build status on TravisCI: [![Build Status](https://travis-ci.com/ibell/pdsim.svg?branch=master)](https://travis-ci.com/ibell/pdsim)

* The documention is available on ReadTheDocs: [![Documentation Status](https://readthedocs.org/projects/pdsim/badge/?version=latest)](http://pdsim.readthedocs.io/en/latest/?badge=latest)

## Install

To setup a conda environment with the requirements, you can do move into the root folder of pdsim and do:
```
conda env create -f conda_environment.yml
```

And then to install PDSim
```
pip -vvv install .
```