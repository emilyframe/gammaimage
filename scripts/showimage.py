# Data Pathway
path = '/Users/eframe/'

import numpy as np
import sys
sys.path.append( path + 'gammaimage/src/')
import functions as func
import tables
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Opening Reconstructed Image Data File
image = path + 'gammaimage/images/image.h5'
outfile = path + 'gammaimage/images/image.gif'

X, Y, Z = np.mgrid[ -49:51:2, -49:51:2, -49:51:2 ]

f = tables.open_file( image, 'r' )
vals = f.root.image.read().reshape( X.shape )
f.close()

x, y = np.mgrid[ -49:51:2, -49:51:2 ]
z = np.arange( -49, 51, 2 )[:]

fig, ax, im = func.make2DMesh( x, y, vals[ :, :, 0 ].T, vmin = min( vals.flatten() ), vmax = max( vals.flatten() ) )
title = plt.title( 'depth = %i mm' %z[ 0 ], fontsize = 20 )

def animate( i ):
    grid = vals[ :, :, i ]
    grid = grid.reshape( x.shape[ 0 ], x.shape[ 1 ] )
    im.set_array( grid.ravel() )
    title.set_text( 'depth = %i mm' %( z[ i ] ) )
    return im

anim = animation.FuncAnimation( fig, animate, frames = len( z ), interval = 400 )
plt.show()

anim.save(outfile)
