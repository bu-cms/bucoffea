from bucoffea.plot.stack_plot import *
import mplhep as hep

# hep styling
plt.style.use(hep.style.CMS)

def move_overflow_to_last_bin(a):
    b = a[:-1]
    b[-1] = b[-1]+a[-1]
    return b

def plot_effparts(acc, region, distribution, distribution_Vmatched, year,  data, mc_all, mc_realV, mc_fakeV, outdir='./output/stack/', integrate=None, ylim=None, xlim=None, rylim=None, tag=None, output_format='pdf', mcscale=1):
    """Creates a data vs MC comparison plot

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    """
    # Rebin
    s = Style()
    h = copy.deepcopy(acc[distribution])
    h_Vmatched = copy.deepcopy(acc[distribution_Vmatched])

    assert(h)
    try:
        newax = hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [240,280,340,400,500,600,1000])
        h = h.rebin(h.axis(newax.name), newax)
        h_Vmatched = h_Vmatched.rebin(h_Vmatched.axis(newax.name), newax)
    except KeyError:
        pass

    all_datasets = list(map(str, h.identifiers('dataset')))
    mapping = {
            'real-V MC' : [ x for x in all_datasets if re.match(mc_realV, x) ],
            'fake-V MC' : [ x for x in all_datasets if re.match(mc_fakeV, x) ],
            }

    mc = re.compile("(fake|real)-V MC")
    for ds in all_datasets:
        if not ( re.match(mc_realV, ds) or re.match(mc_fakeV, ds) ):
            mapping[ds] = [ds];

    h = h.group("dataset", hist.Cat("dataset", "Primary dataset"), mapping)
    h_Vmatched = h_Vmatched.group("dataset", hist.Cat("dataset", "Primary dataset"), mapping)

    # Pick the region we want to look at
    # E.g. cr_2m_j = Di-Muon control region with monojet selection
    h = h.integrate(h.axis('region'),region)
    h_Vmatched = h_Vmatched.integrate(h_Vmatched.axis('region'),region)

    # Plotting
    # Add ratio plot at the bottom if specified (default)
    # Otherwise just plot the histogram
    fig, ax = plt.subplots(1, 1, figsize=(12,10))

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    }
    data_subs_realV_err_opts = {
        #'linestyle':'-',
        'color':'crimson',
        'elinewidth': 1,
    }

    h_data = h.integrate("dataset", data)
    x = h_data.axis('jetpt').centers()
    edges = h_data.axis('jetpt').edges()
    widths = edges[1:] - edges[:-1]
    x = np.append(x, x[-1]+widths[-1])
    edges = np.append(edges, edges[-1]+widths[-1])
    widths = np.append(widths, widths[-1])
    y_data, y_data_e2 = (h_data.values(sumw2=True, overflow='over'))[()]
    y_data    = (y_data)/widths
    y_data_e2 = (y_data_e2)/(widths*widths)
    ax.errorbar(x, y_data, yerr=np.sqrt(y_data_e2), label="data", ls='none', marker='.', markersize=7, color='k', elinewidth=1)
    
    h_mc_fakeV = h.integrate("dataset", re.compile("fake-V MC")).values(sumw2=True, overflow='over')
    y_fakeV, y_fakeV_e2 = h_mc_fakeV[()]
    y_fakeV    = (y_fakeV)/widths
    y_fakeV_e2 = (y_fakeV_e2)/(widths*widths)

    h_mc_realV = h_Vmatched.integrate("dataset", re.compile("real-V MC")).values(sumw2=True, overflow='over')
    y_realV, y_realV_e2 = h_mc_realV[()]
    y_realV    = (y_realV)/widths
    y_realV_e2 = (y_realV_e2)/(widths*widths)
    y_data_subs_realV = y_data - y_realV
    y_data_subs_realV_e2 = y_data_e2 + y_realV_e2

    ax.errorbar(x, y_data_subs_realV, yerr=np.sqrt(y_data_subs_realV_e2), label="data (subtract realV MC)", color="crimson", ls='none', marker='.', markersize=7, elinewidth=1)

    logymin = int(np.log10(max(min(y_data_subs_realV),10e-2))) - 2
    logymax = int(np.log10(min(max(y_data),10e4))) + 1

    # Note the syntax we use to pick the data set
    #hist.plot1d(
    #    h[data],
    #    overlay='dataset',
    #    error_opts=data_err_opts,
    #    ax=ax,
    #    overflow='all',
    #    binwnorm=1)

    hep.histplot(
            [y_fakeV, y_realV],
            edges,
            stack=True,
            label=["fake-V MC", "real-V MC"],
            histtype="fill",
            )


    # Plot MC background samples
    # Here we use a regular expression to match
    # data sets we want
    # if mc!=None:
    #     hist.plot1d(
    #         h[mc],
    #         overlay='dataset',
    #         stack=True,
    #         clear=False,
    #         overflow='all',
    #         ax=ax,
    #         binwnorm=1)

    # Apply correct colors to BG histograms
    handles, labels = ax.get_legend_handles_labels()
    new_labels = []
    for handle, label in zip(handles, labels):
        col = None
        for k, v in colors.items():
            if re.match(k, label):
                col = v
                break
        if col:
            handle.set_color(col)
            handle.set_linestyle('-')
            handle.set_edgecolor('k')

        l = None

        channel = channel_name(region)
        # Pick the proper legend labels for the channel
        if channel == 'VBF':
            legend_labels_to_use = legend_labels['VBF']
        elif channel in ['Monojet', 'Mono-V']:
            legend_labels_to_use = legend_labels['Monojet/Mono-V']

        # Add in the common labels
        legend_labels_to_use.update(legend_labels['Common'])

        for k, v in legend_labels_to_use.items():
            if re.match(k, label):
                l = v
        new_labels.append(l if l else label)

    # Update legend
    try:
        region_name = s.region_names[region]
    except KeyError:
        region_name = region
    ax.legend(title=region_name,ncol=2,handles=handles, labels=new_labels)


    ax.text(1., 0., distribution,
                fontsize=10,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
    fig.text(0., 1., '$\\bf{CMS}$ internal',
                fontsize=14,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

    fig.text(1., 1., f'{channel_name(region)}, {lumi(year):.1f} fb$^{{-1}}$ ({year})',
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
    # Aesthetics
    ax.set_yscale("log")
    ax.set_ylabel('Events / GeV')
    plot_settings=style.plot_settings()
    if region in plot_settings.keys():
        plot_settings=plot_settings[region]
    if distribution in plot_settings.keys():
        plot_settings=plot_settings[distribution]
    if ylim:
        if ylim=="auto":
            width = np.diff([x for x in h.axes() if "dataset" not in x.name][0].edges())
            vmc = h[mc].integrate("dataset").values()[()] / width
            try:
                vdata = h[data].integrate("dataset").values()[()] / width
            except:
                vdata = vmc
            if signal:
                vsig = h[signal].integrate("dataset").values()[()] / width
            else:
                vsig = vmc


            ax.set_ylim(
                0.5*min([np.min(vmc[vmc>0]), np.min(vdata[vdata>0]), np.min(vsig[vsig>0])]),
                1e2*max([np.max(vmc), np.max(vdata), np.min(vsig)]),
            )

        else:
            ax.set_ylim(ylim[0],ylim[1])
    elif 'ylim' in plot_settings.keys():
        ax.set_ylim(plot_settings['ylim'])
    else:
        ax.set_ylim(1e-1,1e6)

    # override ylim
    ax.set_ylim(pow(10, logymin), pow(10, logymax))

    if xlim:
        ax.set_xlim(xlim[0],xlim[1])
    elif 'xlim' in plot_settings.keys():
        ax.set_xlim(plot_settings['xlim'])

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for form in output_format.split(','):
        outpath = pjoin(outdir, f"{region}_{distribution}_{tag + '_' if tag else ''}{year}.{form}")
        fig.savefig(outpath)
        print(f"Saved plot file in {outpath}")
    plt.close('all')
