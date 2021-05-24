# Imported Libraries
import numpy as np
import matplotlib.pyplot as plt

def computeMLEM( sysMat, counts, sens, eps, nIter = 10 ):
    """this function computes iterations of MLEM it returns the image after nIter iterations

    sysMat is the system matrix, it should have shape: (n_measurements, n_pixels)
    It can be either a 2D numpy array, numpy matrix, or scipy sparse
    matrix

    counts is an array of shape (n_measurements) that contains the number
    of observed counts per detector bin

    sens_j is the sensitivity for each image pixel
    if this is None, uniform sensitivity is assumed
    """

    nPix = sysMat.shape[ 1 ]
    lamb = np.ones( nPix )

    iIter = 0
    while iIter < nIter:
        sumKlamb = sysMat.dot( lamb ) + 10e-20
        outSum = ( sysMat * counts[ :, np.newaxis ] ).T.dot( 1 / sumKlamb )
        lamb = ( outSum * lamb ) * ( sens / ( sens ** 2 + max( sens ) ** 2 * eps ** 2 ) )
        iIter += 1

    return lamb

def oversample( a, n ):
    """oversamples bins in an array by repeating the elements n x n times; lengh
    of the array must have an integer square root.

    This is used for oversampling the detector data and the detector efficiency
    data.
    """
    a1 = a.repeat( int( n ), axis = 0 )
    a2 = a1.reshape( int( np.sqrt( len( a ) ) ), int( n * np.sqrt( len( a ) ) ) )
    a3 = a2.repeat( int( n ), axis = 0 )
    return a3

def make2DMesh(x, y, vals, vmin, vmax):
    """generates a 2D mesh plot
    """
    fig, ax = plt.subplots()
    ax.set_xlabel( 'X (mm)', fontsize = 15 )
    ax.set_ylabel( 'Y (mm)', fontsize = 15 )
    ax.tick_params( labelsize = 15 )
    im = ax.pcolormesh( x, y, vals.T, vmin = vmin, vmax = vmax, shading = 'gouraud' )
    cbar = fig.colorbar( im, format = '%.0e' )
    cbar.set_label( label = 'intensity', rotation = 270, fontsize = 15, labelpad = 25 )
    cbar.ax.tick_params( labelsize = 15 )
    return fig, ax, im
