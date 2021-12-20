# a stupid straight-forward fjord model
#
# replace the blocks with an array
# for 1d this works
#
# profiling with
#  python -m cProfile -o prof.log icecoremodel.py 


# blocks[0] is the calving front


import numpy as np
import pylab as plt


class FjordModel(object):

    def __init__(self, nslots):
        self.nslots = nslots
        self.blocks  = np.zeros(nslots, float)
        self.blength = 10
        self.dxrand  = 20
        self.dxcalv  = 10
        self.bias    = 0.3
        self.actblocks = 0


    def spread_blocks(self):
        if self.blocks[0] < 0:
            self.blocks[0] = 0
        for i in range(self.actblocks-1):
            dx = self.blocks[i+1] - self.blocks[i] 
            if dx < self.blength:
                self.blocks[i+1] = self.blocks[i] + self.blength

    def add_blocks(self, newblocks):
        self.blocks[newblocks:] = self.blocks[:-newblocks] + newblocks*self.blength
        self.blocks[:newblocks] = np.arange(newblocks, dtype=float)*self.blength
        self.actblocks += newblocks
        self.actblocks = min(self.actblocks, self.nslots)


nslots  = 500
nblocks = 100
nsteps  = 5000


m = FjordModel(nslots)

m.add_blocks(nslots)
m.spread_blocks()

allpos = np.zeros((nslots, nsteps), float) #+ np.nan

calving = []
# move
for step in range(nsteps):
    dxs = (np.random.random(m.actblocks) - m.bias) * m.dxrand
    for i, dx in enumerate(dxs):
        if i == m.actblocks-1:
            pass
        else:
            if dx > 0:
                dx = max(min(m.blocks[i+1] - m.blocks[i] - m.blength, dx), 0)
            if dx < 0 and (i>0):
                dx = min( -min(m.blocks[i] - m.blocks[i-1] - m.blength, -dx), 0)
        m.blocks[i] += dx

    # calving: add  blocks
    if m.blocks[0] > m.dxcalv:
        print('-----> calving .. ', step)
        m.add_blocks(nblocks)
        # m.spread_blocks()
        calving.append(step)

    m.spread_blocks()
    allpos[-m.actblocks:,step] = m.blocks[:m.actblocks]


front = (np.diff(allpos,axis=0) > m.blength).argmax(axis=0) * m.blength

tt = np.arange(nsteps)
for xs in allpos:
    # plt.plot(xs, tt, 'k.-', lw=0.5)
    plt.plot(xs, tt, 'ks', ms=1)

plt.plot(front, tt, 'red')

stop

tt = np.arange(nsteps)
data = np.diff(allpos, axis=0)[:,::-1].T
#plt.pcolormesh(allpos.T[:,:-1], tt[:, np.newaxis], data, vmin=m.blength, vmax=1.1*m.blength)
times = np.zeros_like(data)
times[:,:] = tt[:,np.newaxis]
#plt.scatter(data, times, c=data, vmin=m.blength, vmax=10*m.blength,)
plt.imshow(data, vmin=m.blength, vmax=3*m.blength)
plt.axis('auto')
