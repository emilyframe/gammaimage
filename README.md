# Introduction
`gammaimage` is a python library for reconstructing gamma-ray images from coded aperture imager data. This library contains scripts for building the coded aperture system response ('sysresponse.py'), performing MLEM image reconstruction ('mlem.py'), and providing a graphical image of the reconstructed image data ('showimage.py').

## Python Dependencies
- numpy
- tables
- sys
- time
- matplotlib.pyplot (if using the 'showimage.py' script)
- matplotlib.animation (if using the 'showimage.py' script)

# Scripts ('gammaimage/scripts')
- 'sysresponse.py': builds the system response and stores the system response as a 2D matrix of i x j elements, where i is the number source voxels in the image space and j is the number of detector pixels. The size of this matrix is currently set to 125000 x 5776. The system response matrix is saved to 'gammaimage/sysresp/sysmat<pose #>.h5'  Furthermore, the system response matrix can be produced at various projection angles (detector poses) for tomographic purposes. If multiple detector poses are given in the 'poses' array, a system response matrix will be saved in a separate file for each pose and labeled accordingly. Currently, only one detector pose is given at 0 degrees.
    - Note: this script reads in data from 'gammaimage/efficiency.h5' and 'gammaimage/farfieldmap.h5'. You will need to download the file 'farfieldmap.h5' from: https://drive.google.com/file/d/11rHENsmg_fgdNAQ7WcaBRLhSCUq6bqV2/view?usp=sharing
- 'mlem.py': performs MLEM image reconstruction using the system response data in 'gammaimage/sysresp/sysmat<pose #>.h5' and the detector data in 'gammaimage/sysresp/data/phantom<pose #>.h5'. This script also allows for tomographic image reconstruction, i.e. incorporating system responses and detector data files from multiple projection angles (detector poses). Currently, however, only one detector pose is provided (at 0 degrees). The final image is saved to 'gammaimage/images/image.h5'.
- 'showimage.py': converts the image data in 'gammaimage/images/image.h5' to a 3D volume and animates 2D slices of this volume over various depths (perpendicular to the imager).
-

## Instructions for Running the Scripts
1. First: build system response matrix. Change data pathway at top of 'gammaimage/scripts/sysresponse.py' script to the appropriate pathway. To run from main directory: 'python scripts/sysresponse.py'
    - Note: it takes roughly 3 minutes to build one system response
1. Second: perform MLEM image reconstruction. Change data pathway at top of 'gammaimage/scripts/mlem.py' script to the appropriate pathway. To run from main directory: 'python scripts/mlem.py'
    - Note: you can adjust the number of iterations by changing the nIter variable. Currently, set to 100.
    - Note: it takes roughly 7 minutes to perform 100 iterations.
1. Third: plot reconstructed data. Change data pathway at top of 'gammaimage/scripts/showimage.py' script to the appropriate pathway. To run from main directory: 'python scripts/showimage.py'
