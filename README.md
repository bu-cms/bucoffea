# BUCoffea

This is HEP analysis code based on the [coffea](https://github.com/CoffeaTeam/coffea) framework.

## Setup

### Python virtual environments
It is generally a good idea to work in a python virtual environment. The virtual environment provides a completely separate python environment to you in which you can easily install all the packages you need and are at the same time separated from other packages installed on your machine. This way, it is easy to exactly define the environment and prevent misconfiguration.

Instructions for how to setup a virtual environment are all over the internet, one example is [here](https://hepdata-lib.readthedocs.io/en/latest/setup.html#sec-setup-virtualenv). 

### Setup on LXPlus / LPC
If you want to run at either LXPlus or the LPC infrastructure, use these instructions to create and activate a working environment:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh
ENVNAME="bucoffeaenv"
python -m venv ${ENVNAME}
source ${ENVNAME}/bin/activate
```

You can leave the environment by typing `deactivate`. You can later reactivate it by sourcing the `activate` script shown above.

### Package installation
This package is installable via pip. For most purposes, you will want to check out a local copy of the repository and then install from that copy after activating your virtual environment.

```
git clone git@github.com:bu-cms/bucoffea.git
python -m pip install -e bucoffea
```

### Running an example
You can run an example analysis to get started:

```
# GRID certificate required for file access!
python ./bucoffea/example/exampleProcessor.py
```

This example will run over two files from two different datasets, produce MET and jet pt histograms, and plot them as PDF files.
Have a look at the code to see how this works under the hood. Also check out the [examples in the Coffea repository](https://github.com/CoffeaTeam/coffea/tree/master/binder), which are much more detailed than this one.


### Running over larger data sets

Currently, HTCondor submission is used to run over large numbers of input files. The submission is implemented in the `execute.py` script.
To get the idea, you run a test job to process some of the SingleMuon data:

```bash
./execute/execute.py -j4  --datasrc 'eos' submit --dataset 'SingleMuon_2017B' --filesperjob 10 --name 'test_submission'
```

This will submit an HTCondor job running on 4 CPUs per node (`-j4`), to run over pre-processed data from my EOS area (`--datasrc eos`). Files related to job submission, as well as the job output can be found in the `submission/*test_submission/` directory. Check it out!

Note that the worker jobs rely on accessing all the executable code from the virtual environment you are using to submit the jobs. Therefore, make sure to have your virtual environment accessible on a shared file system like AFS.

If you are running a larger number of jobs, it's easy to loose track of them. The `monitor.py` script will allow you to track all jobs belonging to a given submission. Usage is easy:

`./execute/monitor.py submission/*test_submission/files`

### XROOTD
Ignore for now.

```
git clone git@github.com:xrootd/xrootd.git
pushd xrootd
mkdir build
pushd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j4
make -j4 install
pip install bindings/python
popd
popd
```
