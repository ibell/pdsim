## 
## Just use docker-compose to spin up this job
##

FROM ubuntu:18.04

RUN apt-get -y -m update && apt-get install -y cmake g++ git zip python-pip 
RUN pip install -U pip

# This ADD forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/ibell/pdsim/git/refs/heads/master pdsim-version.json
RUN git clone --recursive https://github.com/ibell/pdsim

RUN pip install h5py matplotlib cython numpy

ADD https://api.github.com/repos/CoolProp/CoolProp/git/refs/heads/master coolprop-version.json
RUN git clone --recursive https://github.com/CoolProp/CoolProp
WORKDIR /CoolProp/wrappers/Python
RUN python setup.py install

WORKDIR /pdsim
RUN python setup.py install

ENV MPLBACKEND Agg
WORKDIR /pdsim/examples
RUN python simple_example.py
