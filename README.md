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


### Running the analysis

#### Just over a few files

Check out the `./scripts/run_quick_monojet.py` script to see how to quickly run over a few handpicked files, which may be useful for testing. The output will be saved in a file named `monojet_${dataset}.coffea` depending on the name of the dataset you ran over. If you want to look at a few simple histogram plots, you can run a simple plotter with your output:

```
./scripts/print_plots.py monojet_${dataset}.coffea
```

You can check out the code in the script to see how it's done.


#### Using large data sets

Currently, HTCondor submission is used to run over large numbers of input files. The submission is implemented in the `execute.py` script.
To get the idea, you run a test job to process some of the SingleMuon data:

```bash
buexec -j4  --datasrc 'eos' monojet submit --dataset 'SingleMuon_2017B'--async --filesperjob 30 --name "test_submission"
```

This will submit an HTCondor job running on 4 CPUs per node (`-j4`), to run over pre-processed data from my EOS area (`--datasrc eos`). Files related to job submission, as well as the job output can be found in the `submission/test_submission/` directory. Check it out!
The submission script will automatically detect whether you are running on lxplus or at FNAL and will adapt accordingly.

If you are running a larger number of jobs, it's easy to loose track of them. The `bumon` tool will allow you to track all jobs belonging to a given submission. Usage is easy:

`bumon submission/test_submission/`

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
