# Example Coffea Processor

The file `exampleProcessor.py` contains a simple coffea processor example. It will process two job files and then produce some plots under the `tmp_plots` directory. To execute the processor, simply run:

```bash
./exampleProcessor.py
```

## Declaring Histograms

In coffea, one can define histograms using the `coffea.hist` module. In a histogram declaration, one should specify the name, label and the axes of the
histogram as such:

```python
example_hist = hist.Hist("Counts", dataset_axis, met_axis)
```

Axes themselves are also defined with the same module. For binned variables, `hist.Bin` can be used to define an axis with a given range and number of bins: 

```python
met_axis = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 600, 0.25, 1000)
```

For categorical axes (e.g. dataset name), one can use `hist.Cat` instead:

```python
dataset_axis = hist.Cat("dataset", "Primary dataset")
```

## Filling Histograms

Histograms are filled with the `.fill()` method. While calling this method, one should specify values for each of the axes:

```python
output['met'].fill(dataset=dataset,
            met=df["MET_pt"].flatten())
```

## Plotting Histograms

Histograms can be plotted by invoking the `hist.plot1d()` function. To plot a given axis of the histogram, one can integrate over the other axes of the histogram and call the function as follows:

```python
# Get the histogram from the output accumulator
h = output[variable_name]
# Integrate over the dataset name
h = h.integrate('dataset', <dataset_we_want>)

fig, ax = plt.subplots()
hist.plot1d(h, ax=ax)
```

Alternatively, one can plot an axis overlaid by a second axis (e.g. "MET" value for different datasets in the same figure). To do that, simply pass `overlay=<ax_name>` to the function call, where `<ax_name>` is the name of the axis you want to do the overlaying. An example is here:

```python
# Get the histogram from the output accumulator
h = output[variable_name]

# Do not integrate over different datasets

fig, ax = plt.subplots()
# Plot values from different datasets on the same figure
hist.plot1d(h, ax=ax, overlay='dataset')
```

