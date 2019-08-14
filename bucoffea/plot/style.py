from matplotlib import pyplot as plt

def markers(tag):
    if tag =='data':
        ret = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
            'elinewidth': 1,
            'emarker': '_'
        }
    return ret

def matplotlib_rc():
    plt.rc('mathtext',rm='Helvetica')
    plt.rc('mathtext',it='Helvetica')
    plt.rc('mathtext',bf='Helvetica')
    params = {'font.size':14, 'lines.linewidth' : 1}
    plt.rcParams.update(params)