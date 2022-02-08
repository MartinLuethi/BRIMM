# a stupid straight-forward fjord model
#
# replace the blocks in fjordmodel.py with an array
# for 1d this works

# profiling with
# python -m cProfile -o prof.log fjordmodel1d.py


# blocks[0] is the calving front


import numpy as np
import pylab as plt


class FjordModel(object):

    def __init__(self, nslots, param):
        self.nslots  = nslots
        self.param   = param
        self.blength = param.get('blength', 100)
        self.gblength = param.get('gblength', self.blength/2.)
        self.dxrand  = param.get('dxrand', 100)
        self.dxcalv  = param.get('dxcalv', 10)
        self.bias    = param.get('bias', None)
        self.vfront  = param.get('vfront', 0.)

        # initialize blocks
        # np.nan means: no block present
        self.xblocks  = np.zeros(nslots, float) + 1e6 #+ np.nan   
        self.xfront  = 0.    # position of the calving front

    def spread_blocks(self):
        """
        rule
        starting from the calving front, move each block such that
        blocks are never overlapping
        """
        if self.xblocks[0] < self.xfront:       # limit blocks to outside of calving front
            self.xblocks[0] = self.xfront                      
        for i in range(self.nslots-1):
            dx = self.xblocks[i+1] - self.xblocks[i] 
            if dx < self.blength:
                self.xblocks[i+1] = self.xblocks[i] + self.blength

    def calve_blocks(self, calvblocks):
        self.xfront -= calvblocks*self.gblength                          # move calving front back
        self.xblocks[calvblocks:] = self.xblocks[:-calvblocks]              # make room for newly calved blocks, discard the outermost ones
        self.xblocks[:calvblocks] = self.xfront + np.arange(calvblocks)*self.blength      # put the blocks packed

    def add_blocks(self, newblocks):
        """add packed blocks at the calving front at self.xblocks[0:newblocks]"""
        self.xblocks[newblocks:] = self.xblocks[:-newblocks] 
        self.xblocks[:newblocks] = np.arange(newblocks, dtype=float)*self.blength
        # self.actblocks += newblocks
        # self.actblocks = min(self.actblocks, self.nslots)

    def evolution_uncoupled(self, nsteps, nblocks, dxrand=None, bias=None, param=dict()):
        """
        rule
        move a block until it touches the next one
        blocks are never overlapping

        self.xblocks[0] is at the terminus position of the glacier

        bias is an array and pre-defined for every timestep
        """
        if not (dxrand is None):
            self.dxrand = dxrand
        if not (bias is None):
            self.bias = bias
    
        self.add_blocks(nblocks)
        self.spread_blocks()

        self.allpos  = np.zeros((self.nslots, nsteps), float) #+ np.nan
        self.calving_steps = [(0, nsteps)]

        for step in range(nsteps):
            self.xblocks[0] = max(self.xblocks[0], self.xfront)                      # move the first block downstream of the calving front
            self.xblocks[-1] += 1000                                                # move the last block away, otherwise this would be boundary jamming
            self.xblocks[self.xblocks > flength] = flength+5000                      # let the blocks outside of the fjord swim away
            dxs = (np.random.random(self.nslots) - 0.5 + self.bias[step]) * self.dxrand[step]    # moving distance
            if slowzone:
                idx = (slowzone[0] <= self.xblocks) & (self.xblocks <= slowzone[1])
                dxs[idx] *= 0.5
            for i, dx in enumerate(dxs[:-1]):
                if (dx > 0):
                    db = self.xblocks[i+1] - self.blength - self.xblocks[i]
                    self.xblocks[i] += min(db, dx)
                elif (dx < 0) and (i > 0):
                    db = self.xblocks[i] - self.blength - self.xblocks[i-1]
                    self.xblocks[i] += -min(db, -dx)
                elif (dx < 0) and (i == 0):
                    self.xblocks[i] += -min(self.xblocks[i]-self.xfront, -dx)

            # calving criterion: last block moves too far
            if self.xblocks[0] > self.xfront + self.dxcalv:
                print('calving', step)
                self.calve_blocks(nblocks)
                print(len(np.nonzero(self.xblocks < flength+500)[0]))

                # fig, ax = plt.subplots()
                # self.allpos[:,step] = self.xblocks.copy()
                # ax.eventplot(m.allpos.T, colors='r', lineoffsets=1,
                #                     linewidths=3, linelengths=0.8, alpha=0.5)
                # self.spread_blocks()
                # self.allpos[:,step] = self.xblocks.copy()
                # ax.eventplot(m.allpos.T, colors='k', lineoffsets=1,
                #                     linewidths=1.5, linelengths=1)

                # ax.plot(self.xfront, step, 'ms')

                # ax.set_xlim(-100, 2000)
                # plt.show(block=True)

            self.xfront += self.vfront * dt
            self.spread_blocks()

            self.allpos[:,step] = self.xblocks.copy()


# main model run

import os
import xarray

outdir = '../modelruns'
os.makedirs(outdir, exist_ok=True)

nblocks = 20
nsteps  = 20000
nslots  = 600
blength = 50
gblength = 20
flength = 20000     # length of the fjord
dxcalv  = 10
dxrand0 = 20
bias0  = 0.0
slowzone = [10000, 11000]
vfront  = 10    # velocity of calving front in m/day

dt = 0.1       # day
omega = 2*np.pi / 365.

tt = np.arange(nsteps)*dt

# for bias in np.arange(0, 0.41, 0.05):
#     for dxrand in [10, 20, 50, 100]:
for bias0 in [0.1]:
    for dxrand0 in [100]:
        param = dict(nblocks=nblocks, nsteps=nsteps, nslots=nslots,
                     dxcalv=dxcalv, blength=blength,
                     gblength=gblength, flength=flength,
                     dxrand0=dxrand0, bias0=bias0, slowzone=slowzone,
                     dt=dt, vfront=vfront)

        bias   = np.zeros(nsteps) + bias0 
        dxrand = np.zeros(nsteps) + dxrand0
        # sine
        # dxrand = np.zeros(nsteps) + dxrand0 * 0.5*(1 + np.sin(omega*tt))
        # square
        # dxrand = np.zeros(nsteps)+0.2*dxrand0
        # for i in range(10):
        #     dxrand[(2*i)*1000:(2*i+1)*1000] = 1*dxrand0

        m = FjordModel(nslots, param=param)
        m.evolution_uncoupled(nsteps, nblocks, bias=bias, dxrand=dxrand, param=param)

        tt    = np.arange(nsteps)    # times for all timesteps
        front = (np.diff(m.allpos,axis=0) > 1.01*m.blength).argmax(axis=0) * m.blength

        ds = xarray.Dataset(
            {'front':  (['time'], front),
             'blocks': (['slots', 'time'], m.allpos)},
            coords= dict(time=tt, slots=np.arange(nslots)), 
            attrs = param)

        # outfilename = outdir + '/fjordmodel1d_slowzone_{nslots:.0f}_{blength:.0f}_{nblocks:.0f}_{dxcalv:.0f}_{dxrand:.0f}_{bias:.2f}.nc'.format(**param)

        # ds.to_netcdf(outfilename)

fig, ax = plt.subplots()
# box = [[slowzone[0], slowzone[1], slowzone[1], slowzone[0]],
#        [0, 0, 1000, 1000]]
# ax.fill(box[0], box[1], 'm', alpha=0.3)
ax.eventplot(m.allpos.T[::10], colors='k', lineoffsets=1,
                    linewidths=1.5, linelengths=1)
ax.set_xlim(-4000, flength)


stop

# m.xblocks = 5 * blength * np.arange(nslots) + np.random.random(m.nslots)*blength 
# m.spread_blocks()

#m.evolution_uncoupled(nsteps, nblocks)

# plot the results
tt    = np.arange(nsteps)    # times for all timesteps
front = (np.diff(m.allpos,axis=0) > 1.01*m.blength).argmax(axis=0) * m.blength

plt.plot(front, tt, )
#plt.plot(front, tt, 'red')
stop

for t, xs in zip(tt, m.allpos.T):
    plt.plot(xs, t+np.zeros(nslots), 'k.', ms=1)
    # for x in xs:
    #     plt.plot([x,x+95], [t, t], 'k', lw=2)
    #     plt.plot([x+95,x+95], [t-0.2, t+0.2], 'k', lw=0.5)

plt.plot(front, tt, 'red')
plt.xlim(-2000, flength+1000)
stop

for c, n in m.calving_steps:
    t = tt[c]
    plt.plot([0, n*m.blength], [t,t], 'm', lw=3)


plt.figure()
for c0, c1 in zip(m.calving_steps[3:-1][0], m.calving_steps[4:][0]):
    t = tt[c0:c1]
    f = front[c0:c1] 
    plt.plot(f, (t-t[0]), '.-')


stop

tt = np.arange(nsteps)
data = np.diff(m.allpos, axis=0)[:,::-1].T
#plt.pcolormesh(m.allpos.T[:,:-1], tt[:, np.newaxis], data, vmin=m.blength, vmax=1.1*m.blength)
times = np.zeros_like(data)
times[:,:] = tt[:,np.newaxis]
#plt.scatter(data, times, c=data, vmin=m.blength, vmax=10*m.blength,)
plt.imshow(data, vmin=m.blength, vmax=3*m.blength)
plt.axis('auto')
