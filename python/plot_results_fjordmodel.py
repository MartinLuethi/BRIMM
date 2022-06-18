# analyze the model runs
#
# [2022-06-08]
#
# rsync -rvz tinu@mluethi-tux.client.geo.uzh.ch:/home/tinu/projects/fjordmodel/python/fig_modelruns/*.nc .

import numpy as np
import pylab as plt
import xarray
from scipy.stats import linregress

import glob


filenames = sorted(glob.glob('fig_modelruns/*.nc')) 

fjord_length = 50
n_blocks = 300
# fjord_length = 20
# n_blocks = 200

calv_dist = 20.
block_length = 100.

# biass   = np.arange(0.1, 0.111, 0.01)
# dxrands = np.arange(140, 151, 10)
biass   = np.arange(0.01, 0.101, 0.01)
dxrands = np.arange(20, 151, 10)
calvnbs = [20]

pparam = dict()

allres = []
for block_bias0 in biass:
    res = []
    for block_dxrand0 in dxrands:
        # for calv_nblocks in range(10, 51, 10):
        for calv_nblocks in calvnbs:
            param = dict(calv_nblocks=calv_nblocks, fjord_length=fjord_length, n_blocks=n_blocks,
                         calv_dist=calv_dist, block_length=block_length,
                         block_dxrand0=block_dxrand0, block_bias0=block_bias0)
        
            filename = 'fig_modelruns/fjordmodel1d__fjord_length={fjord_length:d}__n_blocks={n_blocks:.0f}__block_length={block_length:.0f}__calv_nblocks={calv_nblocks:.0f}__calv_dist={calv_dist:.0f}__dxrand={block_dxrand0:.0f}__block_bias0={block_bias0:.2f}.nc'.format(**param)

            try:
                ds = xarray.open_dataset(filename)
            except:
                res.append(np.nan)
                continue

            spinup   = 2
            nspinup  = int(spinup*365/ds.dt)

            time     = ds.time.values[nspinup:]
            front    = ds.front.values[nspinup:]
            xfront   = ds.xfront.values[nspinup:]
            # allpos   = ds.blocks.values[:, nspinup:]

            idx   = np.diff(xfront) < -100
            tcalv = time[:-1][idx]

            # average number of days between calving events
            if 0:
                # plt.plot(np.diff(tcalv))
                tcalvmean = np.diff(tcalv)[2:].mean()
                print(tcalvmean)
                res.append(tcalvmean)

            # average propagation speed of the IMWE
            if 0: 
                outname   = 'IMW_speed'
                labeltext = 'IMW propagation speed (km/d)'
                pparam = dict(cmap=plt.cm.RdYlBu_r, vmin=0, vmax=1.6)
                lmin, lmax = 0.3, 0.8

                idx = np.diff(xfront) < -100
                cidx = np.arange(time.shape[0]-1)[idx]

                vs = []
                for i0, i1 in zip(cidx[:-1], cidx[1:]):
                    tt, xx = time[i0+1:i1], front[i0+1:i1]
                    regr = linregress(tt, xx)
                    t0, t1 = time[i0+1], time[i1]
                    x0, x1 = front[i0+1], front[i1]
                    vs.append( -regr.slope )
                    # vs.append( - (x1-x0)/(t1-t0) )

                    # plt.plot([t0,t1], [x0,x1], 'k')
                    # plt.plot(tt, xx, 'k', alpha=0.5)
                    # plt.plot(tt, xx[-1] + (tt-tt[-1])*regr.slope, 'r')

                res.append(np.mean(vs)/1000.)   # in km/day

                # plt.show(block=True)

            # testing the plot -> OK
            if 0: 
                labeltext = 'bias * dxrand'
                res.append(block_bias0 * block_dxrand0)
                continue

            # mean position of blocks < 50 km (fjordlength)
            if 0:
                labeltext = 'mean position of blocks in fjord'
                nbl = allpos[allpos < 50000].mean()
                res.append(nbl)

            # average calving front position 
            if 0:
                labeltext = 'Average position of calving front (m)'
                res.append(xfront.mean())

            # average advance rate of calving front
            if 0:
                outname   = 'front_rate'
                labeltext = 'Average advance rate of calving front (m/d)'
                pparam = dict(cmap=plt.cm.RdBu, vmin=-20, vmax=20)
                lmin, lmax = -7, 7
                try:
                    regr = linregress(time, xfront)
                    r = regr.slope
                except:
                    r = np.nan

                res.append(r)

            # minimum calving front position (20 km fjordlength)
            if 0:
                outname   = 'min_front_pos'
                labeltext = 'Minimum position of calving front (m)'
                res.append(xfront.min())

            # average propagation duration of the IMWE
            if 0:
                outname   = 'prop_dur'
                labeltext = 'Average propagation duration (d)'
                pparam = dict(cmap=plt.cm.RdBu_r, vmin=0, vmax=100)
                lmin, lmax = 10, 30
                idx  = np.diff(xfront) < -100
                cidx = np.arange(time.shape[0]-1)[idx]
                res.append(np.mean(np.diff(time[cidx])))

                # dd = []
                # for i0, i1 in zip(cidx[:-1], cidx[1:]):
                #     t0, t1 = time[i0+1], time[i1]
                #     dd.append((t0-t1))
                # res.append(np.mean(dd))
                # stop

            # average propagation distance of the IMWE
            if 0:
                outname   = 'prop_dist'
                labeltext = 'Average propagation distance (km)'
                pparam = dict(cmap=plt.cm.RdBu, vmin=5, vmax=20)
                lmin, lmax = 12, 15

                idx  = np.diff(xfront) < -100
                cidx = np.arange(time.shape[0]-1)[idx]

                # plt.plot(time, front)

                dd = []
                for i0, i1 in zip(cidx[:-1], cidx[1:]):
                    t0, t1 = time[i0+1], time[i1]
                    x0, x1 = front[i0+1], front[i1]
                    # plt.plot([t0,t1], [x0,x1], 'k')
                    dd.append((x0-x1))
                res.append(np.mean(dd)/1000.)

                # plt.show(block=True)



                # ctime, cfront = time[1:][idx], front[1:][idx]
                # # plt.plot(ctime, cfront, 'ko')
                # plt.plot(time[cidx], front[cidx], 'mo')
                # plt.plot(time[cidx+1], front[cidx+1], 'ro')
                
                # stop

    allres.append(res)

allres = np.array(allres)


def normalize_img(img, val):
    im = (np.clip(img, val, val+1) - val )
    im[im > 0] = 1
    return im

fig, ax = plt.subplots()
pcm = ax.pcolormesh(dxrands, biass, allres, **pparam, )

im = normalize_img(allres, lmin)
cs1 = ax.contour(dxrands, biass, im, levels=[0.5], colors='red', linewidths=3)
im = normalize_img(allres, lmax)
cs2 = ax.contour(dxrands, biass, im, levels=[0.5], colors='red', linewidths=3)
cs  = ax.contour(dxrands, biass, allres, levels = [lmin, lmax], colors='red', linewidths=3, alpha=0)

cbar = fig.colorbar(pcm, label=labeltext)
cs.set_alpha(1)
cbar.add_lines(cs)

# cs2 = clipcontour(allres, lmax)

ax.set_xlabel('dxrand (m)')
ax.set_ylabel('bias')
ax.set_xlim(15,155)
#plt.colorbar(pcm, label='days')

fig.tight_layout()

cnb = int(calvnbs[0])
fig.savefig(f'fig/{outname:s}__nblocks={cnb}.png', dpi=200)
