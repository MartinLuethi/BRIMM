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
        self.param = param
        self.blength = param.get('blength', 100)
        self.dxrand  = param.get('dxrand', 100)
        self.dxcalv  = param.get('dxcalv', 10)
        self.bias    = param.get('bias', 0.)

        # initialize blocks
        # np.nan means: no block present
        self.blocks  = np.zeros(nslots, float) + 1e6 #+ np.nan   

    def spread_blocks(self):
        """
        rule
        starting from the calving front, move each block such that
        blocks are never overlapping
        """
        if self.blocks[0] < 0:
            self.blocks[0] = 0
        for i in range(self.nslots-1):
            dx = self.blocks[i+1] - self.blocks[i] 
            if dx < self.blength:
                self.blocks[i+1] = self.blocks[i] + self.blength

    def add_blocks(self, newblocks):
        """add packed blocks at the calving front at self.blocks[0:newblocks]"""
        self.blocks[newblocks:] = self.blocks[:-newblocks] 
        self.blocks[:newblocks] = np.arange(newblocks, dtype=float)*self.blength
        # self.actblocks += newblocks
        # self.actblocks = min(self.actblocks, self.nslots)


    def evolution_uncoupled(self, nsteps, nblocks, dxrand=None, bias=None, slowzone=None):
        """
        rule
        move a block until it touches the next one
        blocks are never overlapping
        """
        if dxrand:
            self.dxrand = dxrand
        if bias:
            self.bias = bias
    
        self.add_blocks(nblocks)
        self.spread_blocks()

        self.allpos  = np.zeros((self.nslots, nsteps), float) #+ np.nan
        self.calving_steps = [(0, nsteps)]

        for step in range(nsteps):
            self.blocks[0] = max(self.blocks[0], 0)                         # fix the first block
            self.blocks[-1] += 1000
            self.blocks[self.blocks > flength] += 1000                      # let the last blocks swim away
            dxs = (np.random.random(self.nslots) - 0.5 + self.bias) * self.dxrand    # moving distance
            if slowzone:
                idx = (slowzone[0] <= self.blocks) & (self.blocks <= slowzone[1])
                dxs[idx] *= 0.5
            for i, dx in enumerate(dxs[:-1]):
                if (dx > 0):
                    db = self.blocks[i+1] - self.blength - self.blocks[i]
                    self.blocks[i] += min(db, dx)
                elif (dx < 0) and (i > 0):
                    db = self.blocks[i] - self.blength - self.blocks[i-1]
                    self.blocks[i] += -min(db, -dx)
                elif (dx < 0) and (i == 0):
                    self.blocks[i] += -min(self.blocks[i], -dx)

            # calving criterion
            if self.blocks[0] > self.dxcalv:
                self.blocks[nblocks:] = self.blocks[:-nblocks]
                self.blocks[:nblocks] = np.arange(nblocks)*self.blength
                self.spread_blocks()

            self.allpos[:,step] = self.blocks.copy()





    def evolution_uncoupled_old(self, nsteps, nblocks):
        """
        rule
        move a block until it touches the next one
        blocks are never overlapping
        """
        # allpos is used to collect the history for plotting
        self.allpos  = np.zeros((self.nslots, nsteps), float) #+ np.nan
        self.calving_steps = [(0, nsteps)]

        # move
        for step in range(nsteps):
            vs = (np.random.random(self.nslots) - 0.5 + self.bias) * self.dxrand
            self.blocks += vs
            self.spread_blocks()
            stop
            # for i, v in enumerate(vs):
            #     if i == self.nslots-1:
            #         pass
            #     else:
            #         if v > 0:
            #             v = max(min(self.blocks[i+1] - self.blocks[i] - self.blength, v), 0)
            #         if v < 0 and (i>0):
            #             v = min( -min(self.blocks[i] - self.blocks[i-1] - self.blength, -v), 0)
            #     self.blocks[i] += v

            # # calving: add  blocks
            # if self.blocks[0] > self.dxcalv:
            #     print('-----> calving .. ', step)
            #     r = max(np.random.rand(), 0.2)
            #     ncalv = int(r*nblocks)
            #     self.add_blocks(ncalv)
            #     self.calving_steps.append((step, ncalv))
            # self.spread_blocks()
            # # self.allpos[-self.actblocks:,step] = self.blocks[:self.actblocks]
            # self.allpos[:,step] = self.blocks.copy()

        self.calving_steps = np.array(self.calving_steps)

def evolution_overlap(self):
        """
        rule
        blocks may be overlapping
        move them to new position
        then relax: overlapping blocks slide off the blocks they're sitting on
        """

        m.allpos  = np.zeros((m.nslots, nsteps), float) #+ np.nan
        m.calving_steps = [(0, nblocks)]

        m.blocks = np.arange(nslots)*m.blength*1.1 + np.random.random(m.actblocks) * 10
        for step in range(nsteps):
            xs  = m.blocks
            vs  = (np.random.random(m.actblocks) - 0.5 + m.bias) * m.dxrand
            xs += vs

            # plt.plot(np.diff(xs), 'b')
            # plt.plot(xs - np.arange(nslots)*1.1*m.blength)
            # relaxation
            for i in range(10):
                overlap = np.diff(xs) - m.blength
                overlap[overlap > 0] = 0.   # no need to adjust blocks     
                xs[:-1] += 0.3*overlap
                xs[1:] += -0.3*overlap
            # plt.plot(np.diff(xs), 'orange')
            # plt.plot(xs - np.arange(nslots)*1.1*m.blength)

            xs[0] = max(xs[0], 0)
            m.blocks = xs

            # calving: add  blocks
            if m.blocks[0] > m.dxcalv:
                print('-----> calving .. ', step)
                m.add_blocks(nblocks)
                # m.spread_blocks()
                m.calving_steps.append(step)

            m.allpos[-m.actblocks:,step] = m.blocks[:m.actblocks]


# main model run

import os
import xarray

outdir = '../modelruns'
os.makedirs(outdir, exist_ok=True)

nblocks = 20
nsteps  = 20000
nslots  = 300
blength = 50
flength = 20000     # length of the fjord
dxcalv  = 5
dxrand = 100
bias   = 0.35

# for bias in np.arange(0, 0.41, 0.05):
#     for dxrand in [10, 20, 50, 100]:
for bias in [0.1]:
    for dxrand in [100]:

        param=dict(nblocks=nblocks, nsteps=nsteps, nslots=nslots, dxcalv=dxcalv,
                   blength=blength, flength=flength, dxrand=dxrand, bias=bias)

        m = FjordModel(nslots, param=param)
        m.evolution_uncoupled(nsteps, nblocks, dxrand=dxrand, bias=bias, slowzone=[10000, 11000])
        # m.evolution_uncoupled(nsteps, nblocks, dxrand=dxrand, bias=bias)

        tt    = np.arange(nsteps)    # times for all timesteps
        front = (np.diff(m.allpos,axis=0) > 1.01*m.blength).argmax(axis=0) * m.blength

        ds = xarray.Dataset(
            {'front':  (['time'], front),
             'blocks': (['slots', 'time'], m.allpos)},
            coords= dict(time=tt, slots=np.arange(nslots)), 
            attrs = param)

        outfilename = outdir + '/fjordmodel1d_slowzone_{nslots:.0f}_{blength:.0f}_{nblocks:.0f}_{dxcalv:.0f}_{dxrand:.0f}_{bias:.2f}.nc'.format(**param)

        ds.to_netcdf(outfilename)

# fig, ax = plt.subplots()
# ax.eventplot(m.allpos.T[4000:], colors='k', lineoffsets=1,
#                     linewidths=1, linelengths=1)

# m.blocks = 5 * blength * np.arange(nslots) + np.random.random(m.nslots)*blength 
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
plt.xlim(0, flength+1000)
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
