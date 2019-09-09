import numpy as np

def prep(p):

    p.rc('font', family='serif', size=9, serif='STIXGeneral')
    p.rc('mathtext', fontset='stix')
    p.rc('axes', labelsize=9, titlesize=9)
    p.rc('lines', dash_capstyle='round', linewidth=1, color='black', markersize=3)

    p.rc('xtick.major', size=4, width=0.25, pad=2)
    p.rc('xtick.minor', size=2, width=0.25, pad=2)
    p.rc('ytick.major', size=4, width=0.25, pad=2)
    p.rc('ytick.minor', size=2, width=0.25, pad=2)

    p.rc('legend', fontsize=6, handlelength=2, numpoints=1, labelspacing=0.2, borderpad=0.2, fancybox=False, )

    p.rc('savefig', format='pdf')

    p.rc('figure', max_open_warning=0)

    p.autoscale(tight=True)

    p.rc('figure', figsize=(2.75, 2.55))


def post(f, l=False):

    ax = f.gca()

    from matplotlib.ticker import ScalarFormatter

    for pos in ['bottom', 'top', 'right', 'left']:
        ax.spines[pos].set_linewidth(0.5)

    if l:
        l.get_frame().set_linewidth(0.25)
        l.get_frame().set_edgecolor('black')

    figSize = f.get_size_inches()

    ax.set_position([0.25, 0.16, 0.72, 0.72*figSize[0]/figSize[1]])

    if ax.get_yscale() != 'log':
        fmt = ScalarFormatter()
        fmt.set_powerlimits((-2,2))
        ax.yaxis.set_major_formatter(fmt)

    return
