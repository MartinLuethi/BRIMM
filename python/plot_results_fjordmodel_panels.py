# plot the model runs as 4-panel plot
#
# [2022-06-18]

import numpy as np
import pylab as plt
import xarray
from scipy.stats import linregress

# change bias
filenames = [
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=100__block_bias0=0.01.nc',
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=100__block_bias0=0.02.nc',
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=100__block_bias0=0.03.nc',
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=100__block_bias0=0.04.nc',
    ]

# change dxrand
filenames = [
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=50__block_bias0=0.02.nc',
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=70__block_bias0=0.02.nc',
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=100__block_bias0=0.02.nc',
    'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=20__calv_dist=20__dxrand=130__block_bias0=0.02.nc',
    ]

plt.tick_params(labelsize=26)
# fig = plt.figure(figsize=(38, 21), dpi=50)
# axs = fig.subplots(1, 4, sharey=True)

fig = plt.figure(figsize=(38, 21), dpi=50)

for i, filename in enumerate(filenames):
    param = dict(v.split('=') for v in filename[:-3].split('__')[1:])

    # ax = axs[i]
    ax = plt.subplot(1, 4, i + 1)

    ds = xarray.open_dataset(filename)

    spinup   = 4
    nspinup  = int(spinup*365/ds.dt)

    time     = ds.time.values[nspinup:]
    front    = ds.front.values[nspinup:]
    xfront   = ds.xfront.values[nspinup:]
    allpos   = ds.blocks.values[:, nspinup:]

    idx   = np.diff(xfront) < -100
    tcalv = time[:-1][idx]

    idx = np.diff(xfront) < -100
    cidx = np.arange(time.shape[0]-1)[idx]

    ttt = time - spinup*365

    vs = []
    for i0, i1 in zip(cidx[:-1], cidx[1:]):
        tt, xx = ttt[i0+1:i1], front[i0+1:i1]
        regr = linregress(tt, xx)

        # t0, t1 = ttt[i0+1], ttt[i1]
        # x0, x1 = front[i0+1], front[i1]
        # vs.append( -regr.slope )
        # ax.plot([x0,x1], [t0,t1], 'k')

        ax.plot(xx/1000, tt, 'k', alpha=0.5)
        ax.plot((xx[-1] + (tt-tt[-1])*regr.slope)/1000., tt, 'r')
        ax.plot([0, 15], [ttt[i0]]*2, '--', color='orange', lw=2)

    ax.set_xlim(-1, 16)
    ax.grid(alpha=0.4)
    ax.tick_params(labelsize=26)
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))

    if i == 0:
        ax.set_ylabel('Time (days)', fontsize=35)

    ax.set_title('bias={block_bias0} dxrand={dxrand}'.format(**param), fontsize=20)
        
plt.suptitle("Distance from reference (km)", fontsize=35, y=0.08, x=0.507)

#plt.tight_layout()
plt.savefig(f'fig/all.pdf', dpi=200)
