## 
## Just use docker-compose to spin up this job
##

FROM ubuntu:18.04

RUN apt-get -y -m update && apt-get install -y cmake g++ git zip wget bash

# This ADD forces a build (invalidates the cache) if the git repo contents have changed, otherwise leaves it untouched.
# See https://stackoverflow.com/a/39278224
ADD https://api.github.com/repos/ibell/pdsim/git/refs/heads/master pdsim-version.json
RUN git clone --recursive https://github.com/ibell/pdsim

SHELL ["/bin/bash", "-c"]
RUN wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p $HOME/miniconda && \
    export PATH="$HOME/miniconda/bin:$PATH" && \
    conda config --set always_yes yes --set changeps1 no && \
    conda update -q conda && \
    conda create -y -n test-environment python=2.7 numpy matplotlib h5py cython && \
    source activate test-environment && \
    conda install wxpython && \
    python --version && \
    pip install --pre --trusted-host www.coolprop.dreamhosters.com --find-links http://www.coolprop.dreamhosters.com/binaries/Python/ -U --force-reinstall CoolProp

RUN export PATH="$HOME/miniconda/bin:$PATH" && \
    source activate test-environment && \
    python --version && \
    pip install --pre --trusted-host www.coolprop.dreamhosters.com --find-links http://www.coolprop.dreamhosters.com/binaries/Python/ -U --force-reinstall CoolProp
    
WORKDIR /pdsim
RUN export PATH="$HOME/miniconda/bin:$PATH" && \
    source activate test-environment && \
    python --version && \
    python setup.py install

ENV MPLBACKEND Agg
WORKDIR /pdsim/examples
RUN export PATH="$HOME/miniconda/bin:$PATH" && \
    source activate test-environment && \
    python simple_example.py
