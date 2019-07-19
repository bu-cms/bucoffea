# Monojet analysis

This analysis is the first one implemented in this package. It serves as both a tool for the production of physics results, as well as a testing ground for analyzer architecture. 

## Architecture

The analysis logic is implemented in the `monojetProcessor` class, which inherits from the `coffea` processor class. To un-clutter the processor, setup functions, such as the initial creation of histograms, the conversion of the data frame into candidate objects, etc. is outsourced in the `definitions.py` file.


## Configuration

The construction goal for this analyzer is to separate as far as possible the analysis logic from any specific parameter values. The `dynaconf` package is used to provide configuration handling. For now, configuration parameters are stored in a single config file `config.yaml`. The file has a `default` section, where general parameters can be configured for all years (e.g. if the same jet pt cut is applied for all years you want to look at). To accomodate year-dependent choices (e.g. b tag working points), there are additional configuration sections `era2016`, `era2017`, etc, which overwrite the parameter values in the `default` section. Each of the section defines an environment, and the active environment is defined by the `ENV_FOR_DYNACONF` environment variable. The `monojetProcessor` class decides which environment to select based on the data set name it is processing. Currently, the eras are based on the data-taking year, but this setup would also support a finer granularity, for example for time-dependent corrections.

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
