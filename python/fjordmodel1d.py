# 
# BRIMM    Biased Random-Walk Ice Melange Model
#
# a stupid straight-forward fjord model
#
# replace the blocks in fjordmodel.py with an array
# for 1d this works
#
# profiling with
# python -m cProfile -o prof.log fjordmodel1d.py


import numpy as np
import pylab as plt
import xarray
import os


class Plug(object):
    """ a container holding plug information """

    def __init__(self, pos, nsteps):
        self.x0, self.x1 = pos
        self.active = True
        self.lastactive = -1
        self.newpossible = False
        self.bidx = []     # index of the blocks within this plug
        self.bpos = []     # positions of the blocks within this plug
        self.plugged = np.zeros(nsteps)

class FjordModel(object):

    def __init__(self, n_blocks, param):
        """The BRIMM Model

        |xblocks|   is the main data structure holding the x-positions of the blocks
        |n_blocks|    is the maximum number of blocks available
        |xfront|    position of the calving front
        |glacier_velo|   velocity of the calving front

        |fjord_length|   the length of the fjord
        |block_length|   the length (in flow direction) of each block
        |gblock_length|  the length of the block before calving 
        |calv_dist|    distance between front and first block needed to trigger a calving event
        |block_dxrand0|   maximum random block motion at each time step
        |block_bias0|     how much skewed the random walk is

        The random motion at each time step is

        dx = dxrand * (rand(1) - 0.5 + bias) 

        bias = 0    means left/right motion with the same likelyhood, 
        bias = 0.5  means only motion to the right

        """
        self.n_blocks   = n_blocks
        self.param    = param
        self.block_length  = param.get('block_length', 100)
        self.gblock_length = param.get('gblock_length', self.block_length/2.)
        self.block_dxrand0  = param.get('block_dxrand0', 100)
        self.calv_dist   = param.get('calv_dist', 10)
        self.block_bias0    = param.get('block_bias0', 0.)
        self.glacier_velo  = param.get('glacier_velo', 0.)

        self.plugs    = [Plug(pos, nsteps) for pos in param['plugzones']]

        # initialize blocks
        # blocks not in the domain are set to a high coordinate (+1e6)
        # this has the advantage over np.nan, that we have no special case to deal with
        # (x+np.nan would leak into the domain, so one would have to determine valid blocks at every step )
        self.xblocks = np.zeros(n_blocks, float) + 1e6 
        self.xfront  = 0.    # position of the calving front

    def spread_blocks(self):
        """
        Rule:
        starting from the calving front, move each block such that
        blocks are never overlapping
        """
        if self.xblocks[0] < self.xfront:       # limit blocks to outside of calving front
            self.xblocks[0] = self.xfront                      
        for i in range(self.n_blocks-1):
            dx = self.xblocks[i+1] - self.xblocks[i] 
            if dx < self.block_length:
                self.xblocks[i+1] = self.xblocks[i] + self.block_length

    def spread_blocks_overlapping(self):
        """
        TODO: make a working implementation, this is not good yet
        Rule
        starting from the calving front, move each block such that
        """
        if self.xblocks[0] < self.xfront:       # limit blocks to outside of calving front
            self.xblocks[0] = self.xfront                      
        for i in range(self.n_blocks-1):
            dx = self.xblocks[i+1] - self.xblocks[i] 
            if dx < self.block_length:
                self.xblocks[i+1] = self.xblocks[i] + self.block_length

    def calve_blocks(self, calvblocks):
        self.xfront -= calvblocks*self.gblock_length                             # move calving front back
        self.xblocks[calvblocks:] = self.xblocks[:-calvblocks]              # make room for newly calved blocks, discard the outermost ones
        self.xblocks[:calvblocks] = self.xfront + np.arange(calvblocks)*self.block_length      # put the blocks packed

    def calve_blocks_stacked(self, calvblocks):
        # TODO: not working properly
        self.xfront -= calvblocks*self.gblock_length                             # move calving front back
        self.xblocks[calvblocks:] = self.xblocks[:-calvblocks]              # make room for newly calved blocks, discard the outermost ones
        self.xblocks[:calvblocks] = self.xfront + np.abs(np.random.random(calvblocks)-0.5)*np.arange(calvblocks)*self.block_length      # put the blocks packed

    def add_blocks(self, newblocks):
        """add packed blocks at the calving front at self.xblocks[0:newblocks]"""
        self.xblocks[newblocks:] = self.xblocks[:-newblocks] 
        self.xblocks[:newblocks] = np.arange(newblocks, dtype=float)*self.block_length

    def add_blocks_randomly(self, newblocks):
        """add blocks at random positions"""
        pos = np.random.random(newblocks) * self.param['fjord_length']
        pos = sorted(pos)
        self.xblocks[:newblocks] = pos

    def evolution_uncoupled(self, nsteps, calv_nblocks, dxrand=None, bias=None, param=dict()):
        """
        Rule:
        - move a block until it touches the next one
        - blocks are never overlapping

        self.xblocks[0] is closest to the terminus position of the glacier

        |nsteps|  number of time steps
        |calv_nblocks| number of blocks that are released at each calving event
        |bias|    is an array, and pre-defined for every timestep
        |dxrand|  is an array, and pre-defined for every timestep
        """

        if not (dxrand is None):
            self.dxrand = dxrand
        if not (bias is None):
            self.bias = bias
    
        self.add_blocks(calv_nblocks)
        self.spread_blocks()

        self.allpos  = np.zeros((self.n_blocks, nsteps), float) #+ np.nan
        self.xfronts = np.zeros(nsteps)
        self.tcalv   = [(0., 0.)]

        for step in range(nsteps):
            self.xblocks[0] = max(self.xblocks[0], self.xfront)             # move the first block downstream of the calving front
            self.xblocks[-1] += 1000                                        # move the last block away, otherwise this would be boundary jamming
            self.xblocks[self.xblocks > fjord_length] = fjord_length+5000             # let the blocks outside of the fjord swim away
            dxs = (np.random.random(self.n_blocks) - 0.5 + self.bias[step]) * self.dxrand[step]    # random moving distance of blocks
            
            plugidxs = []
            plugposs = []

            # step 1: remember the blocks in the plugareas
            for plug in self.plugs:
                idx = (plug.x0 <= self.xblocks) & (self.xblocks <= plug.x1) 
                plug.bidx = idx
                plug.bpos = self.xblocks[idx]

            # switch on the plugs
            for plug in self.plugs:
                if not plug.active and plug.newpossible:
                    idx = (plug.x0 <= self.xblocks) & (self.xblocks <= plug.x1)
                    if len(np.nonzero(idx)[0]) > 0.7*(plug.x1-plug.x1)/block_length:
                        plug.active = True

            # let the plug work
            for plug in self.plugs:
                if plug.active:
                    idx = (plug.x0 <= self.xblocks) & (self.xblocks <= plug.x1)
                    dxs[idx] *= 0.01
                    plug.plugged[step] = 1
                
            if slowzone:
                idx = (slowzone[0] <= self.xblocks) & (self.xblocks <= slowzone[1])
                dxs[idx] *= 0.5

            for i, dx in enumerate(dxs[:-1]):
                if (dx > 0):
                    db = self.xblocks[i+1] - self.block_length - self.xblocks[i]
                    self.xblocks[i] += min(db, dx)
                elif (dx < 0) and (i > 0):
                    db = self.xblocks[i] - self.block_length - self.xblocks[i-1]
                    self.xblocks[i] += -min(db, -dx)
                elif (dx < 0) and (i == 0):
                    self.xblocks[i] += -min(self.xblocks[i]-self.xfront, -dx)

            # calving criterion: block at calving front moves more than calv_dist away from calving front
            if self.xblocks[0] > self.xfront + self.calv_dist:
                print('calv front')
                self.tcalv.append((step, self.xfront))
                self.calve_blocks(calv_nblocks)
                # print(len(np.nonzero(self.xblocks < fjord_length+500)[0]))     # how many blocks are actually needed in the model
                for plug in self.plugs:
                    plug.newpossible = True

            # second calving criterion: calving happens when within the first 5000m a open lead is more than 1.5 block lengths
            idx = (self.xblocks - self.xfront < 5000)
            if (np.diff(self.xblocks[idx]) > 2.0*self.block_length).any():
                print('calv lead')
                self.tcalv.append((step, self.xfront))
                self.calve_blocks(calv_nblocks)
                # print(len(np.nonzero(self.xblocks < fjord_length+500)[0]))     # how many blocks are actually needed in the model
                for plug in self.plugs:
                    plug.newpossible = True


            # advance the glacier front
            self.xfront += self.glacier_velo * dt

            # make sure blocks are not stock on top of each other
            self.spread_blocks()

            # switch off (and on) the plug
            # step 2: see how much the blocks in the plug area have moved
            for plug in self.plugs:
                dplug = self.xblocks[plug.bidx] - plug.bpos
                if plug.active:
                    print(step, dplug.sum())
                    if (dplug.sum() > 0.9 * ( (plug.x1 - plug.x0)/self.block_length * self.glacier_velo*dt )):
                        plug.active = False
                        plug.lastactive = step
                        plug.newpossible = False
                else:
                    if plug.newpossible:
                        if (dplug.sum() < 0.2 * ( (plug.x1 - plug.x0)/self.block_length * self.glacier_velo*dt )):
                            plug.active = True
                
            self.allpos[:,step] = self.xblocks.copy()
            self.xfronts[step]  = self.xfront


#============================================================
# main model runs
#============================================================

from matplotlib import animation

nyears   = 4
dt       = 0.05      # time step, in units of day
spinup   = 2
nspinup  = int(spinup*365/dt)
nstepyr  = int(365/dt)
nsteps   = int((spinup+nyears)*nstepyr)     # number of time steps

calv_nblocks   = 20.        # number of calved blocks
block_length   = 100.       # length of floating blocks
gblock_length  = 20.        # length of glacier blocks
fjord_length   = 20000.     # length of the fjord
n_blocks       = 200        # maximum number of icebergs that are tracked
calv_dist      = 20.        # distance from terminus of open lead that triggers calving
block_dxrand0  = 20.        # -> set in for loop below
block_bias0    = 0.0       # -> set in for loop below
slowzone = 0 #[7000, 8000]
#plugzones = [[5000, 5500], [7000, 7500]]
plugzones = []
glacier_velo   = 15       # velocity of calving front in m/day

omega = 2*np.pi / 365.

tt = np.arange(nsteps)*dt    # times for all timesteps

outdir = '../modelruns'
os.makedirs(outdir, exist_ok=True)

def make_result_plot(m, outfilename):
    """plot the results like in the paper"""
    # cut off the spinup years
    allpos = m.allpos[:, nspinup:]
    xfronts = m.xfronts[nspinup:]
    tcalv = np.array([t for t, x in m.tcalv if t >= nspinup]) - nspinup

    idx =  (np.diff(allpos, axis=0) > 1.05*m.block_length).argmax(axis=0)    # index of open lead
    xlead = []
    for i, ii in enumerate(idx):
        xlead.append(allpos[ii, i])
    xlead = np.array(xlead)
    xx = xlead - xfronts

    plt.ioff()
    fig, axs = plt.subplots(1, nyears)
    fig.set_figwidth(16)
    fig.set_figheight(10)
    iax = 0
    for n in range(nyears):
        ax = axs[n]
        tn = n * 365
        i0 = n * nstepyr 
        i1 = i0 + nstepyr
        idx = (i0 <= tcalv) & (tcalv <= i1)
        tc  = tcalv[idx]
        for ii0, ii1 in zip(tc[:-1], tc[1:]):
            ax.plot( xx[ii0:ii1]/1000., tt[ii0:ii1]-tn,  'k', alpha=0.7)
            ax.plot( [1, 20], [tt[ii1]-tn]*2,  '--', color='orange')
    for ax in axs:
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 380)
    for ax in axs[1:]:
        ax.set_yticklabels([])

    fig.savefig(outfilename)
    del fig

if 1:
    for block_bias0 in [0.01, 0.02, 0.03, 0.04, 0.05]:
        for block_dxrand0 in [20, 50, 70, 90, 100, 120, 150]:
        for block_dxrand0 in [30, 50, 60, 70, 80, 110, 130, 140]:
            for calv_nblocks in [10, 20, 30, 40, 50]:
                print('---- running --- :', block_bias0, block_dxrand0, calv_nblocks)
                param = dict(calv_nblocks=calv_nblocks, nsteps=nsteps, n_blocks=n_blocks,
                             calv_dist=calv_dist, block_length=block_length,
                             gblock_length=gblock_length, fjord_length=fjord_length,
                             block_dxrand0=block_dxrand0, block_bias0=block_bias0, 
                             slowzone=slowzone, plugzones=plugzones,
                             dt=dt, glacier_velo=glacier_velo)

                outfilename = 'fig_modelruns/fjordmodel1d__n_blocks={n_blocks:.0f}_block_length={block_length:.0f}_calv_nblocks={calv_nblocks:.0f}_calv_dist={calv_dist:.0f}_dxrand={block_dxrand0:.0f}_block_bias0={block_bias0:.2f}.png'.format(**param)
                
                print(outfilename)
                bias   = np.zeros(nsteps) + block_bias0 
                dxrand = np.zeros(nsteps) + block_dxrand0

                ## sine forcing
                # dxrand = np.zeros(nsteps) + block_dxrand0 * 0.5*(1 + np.sin(omega*tt))
                ## square forcing
                # dxrand = np.zeros(nsteps)+0.2*block_dxrand0
                # for i in range(10):
                #     dxrand[(2*i)*1000:(2*i+1)*1000] = 1*block_dxrand0

                # initialize the calving model
                m = FjordModel(n_blocks, param=param)
                m.add_blocks_randomly(n_blocks)
                m.spread_blocks()

                m.evolution_uncoupled(nsteps, calv_nblocks, bias=bias, dxrand=dxrand, param=param)

                front = (np.diff(m.allpos,axis=0) > 1.01*m.block_length).argmax(axis=0) * m.block_length

                ds = xarray.Dataset(
                    {'front':  (['time'], front),
                     'xfront':  (['time'], m.xfronts),
                     'blocks': (['slots', 'time'], m.allpos)},
                    coords= dict(time=tt, slots=np.arange(n_blocks)), 
                    attrs = param)

                make_result_plot(m, outfilename)
                ds.to_netcdf(outfilename.replace('.png','.nc'))

                stop
stoooop

# iceberg animation 
if 0:
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    step0 = 0
    dstep = 5
    x0, x1 = -2500, 11000
    ax = plt.axes(xlim=(x0, x1), ylim=(0., 2))
    xs = [x0, x1, x1, x0]
    ys = [0.5, 0.5, 1.5, 1.5]
    ax.fill(xs, ys, color='b', alpha = 0.6)
    eplot, = ax.eventplot(m.allpos.T[0], colors='white', lineoffsets=1,
                    linewidths=3, linelengths=1)
    xf = m.allpos[0,0]
    xs = [x0, xf, xf, x0]
    glacier, = ax.fill(xs, ys, 'lightblue', alpha=0.7)
    
    pplugs = []
    for plug in m.plugs:
        xp0, xp1 = plug.x0, plug.x1
        xp = [xp0, xp1, xp1, xp0]
        yp = [0.3, 0.3, 1.7, 1.7]
        pplugs.append( ax.fill(xp, yp, 'm', alpha=0.7, zorder=-1)[0] )

    # initialization function: plot the background of each frame
    def init():
        eplot.set_positions([])
        glacier.set_xy(np.vstack((xs, ys)).T)
        # for pplug in pplugs:
        #     pplug.set_xy(np.vstack((xp, yp)).T)

    # animation function.  This is called sequentially
    def animate(step):
        sstep = step0+dstep*step
        eplot.set_positions(m.allpos[:,sstep])
        xf = m.allpos[0,sstep]
        xs = [x0, xf, xf, x0]
        ys = [0.5, 0.5, 1.5, 1.5]
        a  = np.vstack((xs, ys)).T
        glacier.set_xy(a)

        for plug, pplug in zip(m.plugs, pplugs):
            xp0, xp1 = plug.x0, plug.x1
            xp = [xp0, xp1, xp1, xp0]
            yp = [0.3, 0.3, 1.7, 1.7]
            if plug.plugged[sstep]:
                pplug.set_xy(np.vstack((xp, yp)).T)
            else:
                pplug.set_xy(np.vstack(([0,0,0,0], [0,0,0,0])).T)

    # call the animator.  blit=True means only re-draw the parts that have changed.
    # anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                                frames=int((nsteps-step0)/dstep), interval=20, blit=True)
    anim = animation.FuncAnimation(fig, animate, init_func=init, repeat=False,
                                   frames=int(nsteps/dstep), interval=20) #, blit=True)

    # anim.save('fjordmodel_plug.avi', fps=30)  #, extra_args=['-vcodec', 'libx264'])

    plt.show()
    stop



stop

for t, xs in zip(tt, m.allpos.T):
    plt.plot(xs, t+np.zeros(n_blocks), 'k.', ms=1)
    # for x in xs:
    #     plt.plot([x,x+95], [t, t], 'k', lw=2)
    #     plt.plot([x+95,x+95], [t-0.2, t+0.2], 'k', lw=0.5)

plt.plot(front, tt, 'red')
plt.xlim(-2000, fjord_length+1000)
stop

for c, n in m.calving_steps:
    t = tt[c]
    plt.plot([0, n*m.block_length], [t,t], 'm', lw=3)


plt.figure()
for c0, c1 in zip(m.calving_steps[3:-1][0], m.calving_steps[4:][0]):
    t = tt[c0:c1]
    f = front[c0:c1] 
    plt.plot(f, (t-t[0]), '.-')


stop

tt = np.arange(nsteps)
data = np.diff(m.allpos, axis=0)[:,::-1].T
#plt.pcolormesh(m.allpos.T[:,:-1], tt[:, np.newaxis], data, vmin=m.block_length, vmax=1.1*m.block_length)
times = np.zeros_like(data)
times[:,:] = tt[:,np.newaxis]
#plt.scatter(data, times, c=data, vmin=m.block_length, vmax=10*m.block_length,)
plt.imshow(data, vmin=m.block_length, vmax=3*m.block_length)
plt.axis('auto')
