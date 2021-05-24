# Data Pathway
path = '/Users/eframe/'

# Libraries
import numpy as np
import tables
import sys
sys.path.append( path + 'gammaimage/src/' )
import functions as func
import time

# Input Files
detectorImages = np.array( [ path + 'gammaimage/data/phantom0.h5' ] )   # integer indicates pose in degrees, e.g. 0 degrees
systemMatrices = np.array( [ path + 'gammaimage/sysresp/sysmatrix0.h5' ] ) # should be same length as detectorImages

# Initializing Source Voxels
sourceX, sourceY, sourceZ = np.mgrid[ -49:51:2, -49:51:2, -49:51:2 ]
sourcePixels = np.array( [ sourceX.flatten(), sourceY.flatten(), sourceZ.flatten() ] ).T

# Initializing Detector Pixels
detectorX, detectorY = np.mgrid[ -37.5:38.5:1, -37.5:38.5:1 ]
detectorPixels = np.array( [ detectorX.flatten(), detectorY.flatten() ] ).T

# Loading Detector Data
detectorData = np.zeros( ( len( detectorImages ), len( detectorPixels ) ) )
for i in np.arange( len( detectorImages ) ):
    detectorData[ i ] = tables.open_file( detectorImages[ i ], 'r' ).root.counts.read()

# Loading System Matrix Data
sensitivity = np.zeros( len( sourcePixels ), dtype = 'float32' )
systemMatrix = { str( 0 ): np.zeros( ( len( sourcePixels ), len( detectorPixels ) ), dtype = 'float32' ) }
for i in np.arange( len( systemMatrices) ):
    matrix = tables.open_file( systemMatrices[ i ], 'r').root.matrix.read()
    sensitivity = sensitivity + matrix.sum( axis = 1 )
    systemMatrix[ str( i ) ] = matrix

# Performing MLEM Image Reconstruction
eps = 10e-2                                 # sensitivity weighting factor
nIter = 100                                 # number of iterations
lamb = np.ones( len( sourcePixels ) )
tme = time.time()
for i in np.arange( nIter ):
    print( 'iteration: %i of %i' %( i + 1, nIter ) )
    ratio = np.zeros( ( len( sourcePixels ), len( detectorPixels )  ), dtype = 'float32' )
    for p in np.arange( len( detectorImages ) ):
        sysMat = systemMatrix[ str( p ) ]
        data = detectorData[ p ]
        projExpected = np.dot( sysMat.T, lamb )
        ratio = ratio + sysMat * np.divide( data, projExpected, out = np.zeros_like( data ), where = projExpected != 0 )
    lamb = lamb * ratio.sum( axis = 1 ) * ( sensitivity / ( sensitivity ** 2 + max( sensitivity ) ** 2 * eps ** 2 ) )
print( time.time() - tme )
f = tables.open_file( path + 'gammaimage/images/image.h5', 'w')
f.create_array('/', 'image', lamb)
f.close()
