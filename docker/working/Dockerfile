## 
## Just use docker-compose to spin up this job
##

FROM continuumio/miniconda

RUN apt-get -y -m update && apt-get install -y cmake g++ gfortran git zip gtk2.0

# This ADD forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/ibell/pdsim/git/refs/heads/master pdsim-version.json
RUN git clone --recursive https://github.com/ibell/pdsim

ADD https://api.github.com/repos/CoolProp/CoolProp/git/refs/heads/master coolprop-version.json
RUN git clone --recursive https://github.com/CoolProp/CoolProp

RUN conda install -y wxpython numpy h5py matplotlib cython

WORKDIR /CoolProp/wrappers/Python
RUN python setup.py install

WORKDIR /pdsim
RUN python setup.py install

WORKDIR /pdsim/examples
RUN python simple_example.py
