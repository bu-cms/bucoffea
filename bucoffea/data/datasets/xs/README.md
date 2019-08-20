### Sample cross sections

Sample cross sections are accumulated in two steps.

#### GenXSecAnalyzer

The `GenXSecAnalyzer` is used to determine the cross sections as defined by the generator used to create the MC sample. The `run_auto_genxsec.sh` script does all the heavy lifting of running the analyzer. It will automatically run over all datasets in the `datasets_201*.txt` lists in the `datasets` folder. Please note that you will need to a) **not** source the usual bucoffea python3 environment and b) source a recent CMSSW version for this to work. The script will automatically parse the output of the analyzer and dump it into the `xs.txt` file. Note: The script will only analyze datasets **that are not yet included in xs.txt**, so that after adding one new dataset, you will not need to run over all datasets again. However, this also means that if you want to re-run over a dataset that has already been run over in the past (e.g. to run over a larger number of events), you will need to **remove the datasets from xs.txt first**.

#### Conversion to YAML

To make further steps easier, the raw `xs.txt` list is converted into a YAML file using `xs_to_yaml.py`. The output in `xs.yml` is a nested dictionary, with top-level keys being the short names of the dataset (i.e. the dataset names **we** use rather than the one **CMS** uses.). For each dataset, multiple cross-section values can be stored under arbitrary keys, with the exception of the `gen` key, which is reserved for the result of the GenXSecAnalyzer. Example:

```yaml
# (...)
ZJetsToNuNu_HT-100To200-mg_2017:
  gen: 306.2
ZJetsToNuNu_HT-100To200-mg_2018:
  gen: 303.6
ZJetsToNuNu_HT-1200To2500-mg_2017:
  gen: 0.3421
# (...)
```

The conversion script will automatically fill the `gen` cross sections for each dataset. If we have additional higher-order cross-sections (e.g. NNLO for `ttbar`), we can add them to the yaml file manually. Note that the conversion script will not overwrite our manually added keys. It will only ever touch the `gen` keys.


#### Reading of the YAML cross sections

Technically, this part does not happen here, but in `plot/util.py`, but it's closely related, so it's also documented here. The logic behind reading the cross-sections is the following: If the user does not specifically say which cross section he wants for a given dataset, we take the one with the highest order we can find. For example:

```yaml
# (...)
ZJetsToNuNu_HT-100To200-mg_2017:
  gen: 306.2
  nlo: 350.4
  nnlo: 363.8
# (...)
```

Here, we want to use the `nnlo` cross section by default because it is the highest perturbative order. Note that no matter how the sample was generated, we always consider the `gen` cross section to be the lowest order, i.e. if another cross section is present, we will always prefer that one over the `gen` value.

In specific cases, it one may want to explicitly choose the cross section value. For this case, we introduce the `use` key, which define the cross section we should use. Example:

```yaml
# (...)
ZJetsToNuNu_HT-100To200-mg_2017:
  use: gen
  gen: 306.2
  nlo: 350.4
  nnlo: 363.8
# (...)
```

In this case, we detect that the `use` key is present and check it to find out which cross section to use. Since it says `gen`, we will use that one.

