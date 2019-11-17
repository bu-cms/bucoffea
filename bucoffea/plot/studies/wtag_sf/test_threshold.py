from bucoffea.plot.stack_plot import *
from klepto.archives import dir_archive


inpath = "../../input/merged"
region = 'cr_2m_inclusive_v'
year=2017
mc = re.compile(f'(VDY.*HT.*|QCD.*|W.*HT.*|GJets.*HT.*|ZJetsToNuNu.*){year}')

# load the archive
acc = dir_archive(
    inpath,
    serialized=True,
    compression=0,
    memsize=1e3,
    )
acc.load('sumw')
acc.load('sumw_pileup')
acc.load('nevents')


for distribution in ['ak8_wvsqcd0','ak8_wvsqcdmd0']:
    print(f'Variavle: {distribution}')
    acc.load(distribution)
    acc[distribution] = merge_extensions(acc[distribution], acc, reweight_pu=not ('nopu' in distribution))
    scale_xs_lumi(acc[distribution])
    acc[distribution] = merge_datasets(acc[distribution])
    
    h = acc[distribution]
    h = h.integrate('dataset', mc)

    edges = h.axis('tagger').edges()
    values = h.values()[(region,)]
    ntotal = sum(values)

    print('score cut\t mistag rate')
    for ibin in range(len(edges)):
        cut=edges[ibin]
        ncut = sum(values[ibin:])
        rate = 1.0*ncut/ntotal
        print(f'{cut}\t {rate}')

