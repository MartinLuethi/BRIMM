# analyze the model runs
#
# get data with
# rsync -rvz tinu@mluethi-tux.client.geo.uzh.ch:/home/tinu/projects/fjordmodel/modelruns/*.nc .

import numpy as np
import pylab as plt
import xarray
import glob

filenames = sorted(glob.glob('../modelruns/*_0.05.nc')) 

c_mean, c_max = [], []
dxrands = []
for filename in filenames:
    print(filename)
    ds = xarray.open_dataset(filename)
    time  = ds.time.values
    param = ds.attrs

    front   = ds.front.values
    dfront  = np.diff(front)
    calvidx = (dfront > 0.7*dfront.max())
    calving = time[1:][calvidx]

    # print(np.diff(calving).mean())
    c_mean.append(np.diff(calving).mean())
    c_max.append(np.diff(calving).max())
    dxrands.append(param['dxrand'])

    if 1:
        fig, ax = plt.subplots()
        # box = [[slowzone_start, slowzone_end, slowzone_end, slowzone_start],
        #        [0, 0, 1000, 1000]]
        # ax.fill(box[0], box[1], 'm', alpha=0.3)
        ax.eventplot(ds.blocks.T[-5000::10], colors='k', lineoffsets=1,
                            linewidths=1.5, linelengths=1)
        ax.set_xlim(0, flength)
        stop

    if 0:
        # for i in range(len(calving)-1):
        for i in range(-9,-1):
            i0, i1 = calving[i], calving[i+1]
            dfront = front[i0:i1]-front[i1]
            plt.plot(dfront - dfront[-1], time[i0:i1]-time[i0])
            plt.xlabel('distance from front (m)')
            plt.ylabel('Time (some units)')
        stop

c_mean = np.array(c_mean)
dxrands = np.array(dxrands)
idx = dxrands.argsort()

if 1:
    
    # plt.plot(param['dxrand'], dfront[calvidx].mean(), 'o')
    plt.plot(dxrands[idx], c_mean[idx], '-o')
    #plt.plot(param['dxrand'], c_max, 's')

    plt.xlabel('dxrand')
    plt.ylabel('c_mean')

if 0:
    plt.plot(np.diff(calving))

    # plt.plot(front)
    # plt.plot(calving, np.zeros_like(calving), 'ro')

if 0:
    fig, ax = plt.subplots()
    # box = [[slowzone_start, slowzone_end, slowzone_end, slowzone_start],
    #        [0, 0, 1000, 1000]]
    # ax.fill(box[0], box[1], 'm', alpha=0.3)
    ax.eventplot(m.allpos.T[-1000:], colors='k', lineoffsets=1,
                        linewidths=1.5, linelengths=1)
    ax.set_xlim(0, flength)

