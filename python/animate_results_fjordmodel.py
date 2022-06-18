# iceberg animation 
#
# plot the data of the netcdf files as animation

import numpy as np
import pylab as plt
import xarray
from matplotlib import animation

nyears   = 4
dt       = 0.05      # time step, in units of day
spinup   = 2
nspinup  = int(spinup*365/dt)
nstepyr  = int(365/dt)
nsteps   = int((spinup+nyears)*nstepyr)     # number of time steps


filename = 'fig_modelruns/fjordmodel1d__fjord_length=50__n_blocks=300__block_length=100__calv_nblocks=10__calv_dist=20__dxrand=100__block_bias0=0.01.nc'

ds = xarray.open_dataset(filename)
spinup   = 2
nspinup  = int(spinup*365/ds.dt)

time     = ds.time.values[nspinup:]
front    = ds.front.values[nspinup:]
xfront   = ds.xfront.values[nspinup:]

allpos   = ds.blocks.values[:,nspinup:]

if 1:
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    step0 = 0
    dstep = 20
    x0, x1 = 5000, 51000
    ax = plt.axes(xlim=(x0, x1), ylim=(0., 2))
    xs = [x0, x1, x1, x0]
    ys = [0.5, 0.5, 1.5, 1.5]
    ax.fill(xs, ys, color='b', alpha = 0.6)
    eplot, = ax.eventplot(allpos.T[0], colors='white', lineoffsets=1,
                    linewidths=3, linelengths=1)
    xf = allpos[0,0]
    xs = [x0, xf, xf, x0]
    glacier, = ax.fill(xs, ys, 'lightblue', alpha=0.7)
    
    xm = xfront[0] + front[0]
    mfront, = ax.plot([xm, xm], [0.4, 1.6], 'r', lw=3)

    # initialization function: plot the background of each frame
    def init():
        eplot.set_positions([])
        glacier.set_xy(np.vstack((xs, ys)).T)
        mfront.set_xdata([xm, xm])
        # for pplug in pplugs:
        #     pplug.set_xy(np.vstack((xp, yp)).T)

    # animation function.  This is called sequentially
    def animate(step):
        sstep = step0+dstep*step
        eplot.set_positions(allpos[:,sstep])
        xf = allpos[0,sstep]
        xs = [x0, xf, xf, x0]
        ys = [0.5, 0.5, 1.5, 1.5]
        a  = np.vstack((xs, ys)).T
        glacier.set_xy(a)
        xm = xfront[sstep] + front[sstep]
        mfront.set_xdata([xm, xm])



    # call the animator.  blit=True means only re-draw the parts that have changed.
    # anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                                frames=int((nsteps-step0)/dstep), interval=20, blit=True)
    anim = animation.FuncAnimation(fig, animate, init_func=init, repeat=False,
                                   frames=int(nsteps/dstep), interval=20) #, blit=True)

    # anim.save('fjordmodel_plug.avi', fps=30)  #, extra_args=['-vcodec', 'libx264'])

    plt.show()
    stop
