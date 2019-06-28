# Monojet analysis


## Configuration

The construction goal for this analyzer is to separate as far as possible the analysis logic from any specific parameter values. The `dynaconf` package is used to provide configuration handling. For now, configuration parameters are stored in a single config file `config.yaml`. The file has a `default` section, where general parameters can be configured for all years (e.g. if the same jet pt cut is applied for all years you want to look at. To accomodate year-dependent choice (e.g. b tag working points), there are additional configuration sections `era2016`, `era2017`, etc, which overwrite the parameter values in the `default` section. Each of the section defines an environment, and the active environment is defined by the `ENV_FOR_DYNACONF` environment variable. Currently, this is just set manually, and we should implement an automatization based on the input data set.

### Parameter access

It's easiest to understand how to use the configuration from a simple example. Let's say our configuration file has this part in it:

```yaml
era2016:
  selection:
    signal:
      recoil: 250    # min
```

The logic here is that all cuts for event selection go under the "selection" header. Cuts for the signal selection go under "signal", and the cut on the recoil variable is 250 GeV. We could also specify minimal and maximal cuts with a further nested layer, but this suffices for now.

To read this configuration in python, we do something along the lines of:

```python
os.environ["ENV_FOR_DYNACONF"] = "era2016"
os.environ["SETTINGS_FILE_FOR_DYNACONF"] = os.path.abspath("config.yaml")
from dynaconf import settings as cfg

# Configuration parameters
print(cfg.SELECTION.SIGNAL.RECOIL)
# --> 250
```

The parameters are trivially accessed by using the configuration dictionary keys as attributed to the configuration object. They are given in all caps in order to visually separate them from "normal" python code.

<!-- ## Speed benchmarking

Running the monojet analysis without any systematic variations or external weight lookup on the a TTBar sample with 142M events over XRootD takes X minutes.

`/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19_ext1-v1/NANOAODSIM`
 -->

## 