# a stupid straight-forward fjord model
#
# just blocks stacked on each other
#
# [2021-11-03]
# This is a stripped-down 1D version. Nicely enough, we can use a list
# to store the blocks, which makes handling much simpler
#
# [2021-11-04] this is a pure 2D version, maybe it is faster

# profiling with
#  python -m cProfile -o prof.log icecoremodel.py 


import numpy as np
import pylab as plt

rhog    = 917 * 9.81 * 1e-6    # in MPa


class Block(object):
    """
    k            index in the 1D grid, used to identifiy the block
    lbh          (laenge, breite, hoehe)
    stress:      sig
    strain:      eps
    position:    pos
    rel. density D
    """

    def __init__(self, nr, pos=0, h=1):
        self.nr = nr
        self.bout = None
        self.bin  = None
        self.pos = pos

    def __repr__(self):
        return 'Block({} - {:.1f})'.format(self.nr, self.pos)


blocks  = [Block(0, 0)]
blength = 10
dxrand  = 100
dxcalv  = 10
nblocks = 50
nsteps  = 5000
bias    = 0.2


def spread_blocks(blocks):
    for b in reversed(blocks):
        if b.bin:
            if b.pos - b.bin.pos < blength:
                b.pos = b.bin.pos + blength

def add_blocks(blocks, nblocks):
    nr = blocks[-1].nr
    for i in range(nblocks):
        b = Block(i+nr+1)
        blocks.append(b)

    for b1, b2 in zip(blocks[:-1], blocks[1:]):
        b2.bout = b1
        b1.bin  = b2

add_blocks(blocks, 200)
spread_blocks(blocks)

calving = blocks[0]

allpos = []
calving = []
# move
for step in range(nsteps):
    dxs = (np.random.random(len(blocks)) - bias) * dxrand
    for b, dx in zip(blocks, dxs):
        if dx > 0 and b.bout:
            dx = max(min(b.bout.pos - b.pos - blength, dx), 0)
        elif dx < 0 and b.bin:
            # print('\n',dx)
            dx = min( - min(b.pos - b.bin.pos - blength, -dx), 0)
            # print(dx)
            # print()
        b.pos += dx
        if b.pos < 0:
            b.pos = 0

    # calving: add 20 blocks
    if blocks[-1].pos > dxcalv:
        calving.append(step)
        add_blocks(blocks, nblocks)
        spread_blocks(blocks)

    # # delete blocks that are too far away
    # for i, b in enumerate(reversed(blocks)):
    #     if b.bin:
    #         if (b.pos - b.bin.pos) > 20:
    #             print(b, b.bin)
    #             blocks = blocks[:i]

    allpos.append([b.pos for b in blocks])

    
allallpos = []
a = allpos[0]
aa = [a]
ll = len(a)

for a in allpos:
    if len(a) == ll:
        aa.append(a)
    else:
        allallpos.append(aa)
        aa = [a]
        ll = len(a)

maxy = np.diff(calving).max()

for aa in allallpos:
    aa = np.array(aa)
    xs = np.arange(aa.shape[0])
    for row in aa.T:
        plt.plot(row, xs, 'k-', lw=0.5)
    plt.xlim(0,400*blength)
    plt.ylim(0,maxy+5)
    plt.show(block=True)
stop        

# for i, a in enumerate(allpos):
#     plt.plot(a, [i]*len(a), 'k.')
stop
allpos = np.array(allpos)

xs = np.arange(allpos.shape[0])
for row in allpos.T:
    plt.plot(row, xs, 'k.-', lw=0.5)

