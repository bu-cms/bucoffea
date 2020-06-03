# Dataset management

The overall workflow used in this framework can be summarized as:

1. Make list of interesting data sets.
2. Use CRAB to run nanoaod-tools over the list of data sets to create modified nanoaod "ntuples".
3. Use the coffea processors to analyze the ntuples.

The files in this folder pertain to point 1. above. Below, the details of this step are explained.

## Making data set lists

CMS data sets are generally formatted as 

```
/NAME/CONDITIONS/DATATIER
```

IN this framework, the DATATIER is always NANOAOD for data, and NANOAODSIM for MC. The define a list of interesting datasets for a given NanoAOD campaign, a list of NAME values is defined. These are located in `dataset_names_*txt`, and look like this:

```
# W NLO
W*JetsToLNu_LHEWpT_*_Tune*_13TeV-amcnloFXFX-pythia8
WJetsToLNu_*J_Tune*_13TeV-amcatnloFXFX-pythia8
```

Asterisks represent wildcards. This list of NAME values is converted into properly formatted datasets using the `expand_datasets.sh` script. An `expand` function is defined to be usable e.g. like this:

```
INFILE="dataset_names_mc.txt"
CONDITIONS="RunIISummer16*01June*"
expand "${INFILE}" "{CONDITIONS}" | tee datasets_2016.txt
```

For each dataset NAME stub given in the input files, a dataset name is constructed by combining the NAME and CONDITIONS strings, finding all matching data sets in DAS and then printing the results to the terminal, where they are redirected to a file. Note that the condition string shown here contains the identifier `1June`, which corresponds to NanoAOD v5.

### Adding a new data set

In order to add a new data set, simply add its name to the approriate `dataset_names_*txt` file and rerun the expansion script.

### Updating the NanoAOD campaign

To update the campaign, simply change the conditions strings in the expansion script. NanoAOD campaigns have characteristic identifiers like the `1June` one shown above, so you simply need to find out the one for your target campaign.