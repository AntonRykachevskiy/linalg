import numpy as np
import skimage.filters as filt
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import sobel
from skimage import transform
from skimage import data, io, filters
from skimage import img_as_ubyte

from TILT import *

if __name__ == '__main__':
    check = io.imread('cm_1.jpg')
    print np.amax(check)
    n = min(check.shape[0], check.shape[1])
    check = check[:n, :n]



    init_points = np.asarray([[0, n], [0, n]])

    print "shape",  check.shape

    plt.imshow(check)
    plt.show()
    checkie = np.zeros((n,n,3))
    checkie[:,:,0] = check
    checkie[:,:,1] = check
    checkie[:,:,2] = check
    Ds, Dotau, A, E = TILT(checkie, 'euclidean', init_points)


    plt.imshow(A)
    #plt.imshow(E)
    plt.show()
