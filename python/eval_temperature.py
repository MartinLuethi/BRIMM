# try to figure out a thermodynamic temperture from mean particle mobility / speed

import numpy as np
import pylab as plt
import xarray

ds = xarray.open_dataset('fig_modelruns/fjordmodel1d__n_blocks=200_block_length=100_calv_nblocks=30_calv_dist=20_dxrand=100_block_bias0=0.01.nc')



spinup   = 2
nspinup  = int(spinup*365/ds.dt)

allpos = ds.blocks.values[:, nspinup:]

plt.plot(allpos[0])

