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

check = io.imread('tiny_check.jpg')
print np.amax(check)
init_points = np.asarray([[0, 50], [0, 50]])

checkie = np.zeros((50,50,3))
checkie[:,:,0] = check
checkie[:,:,1] = check
checkie[:,:,2] = check
Ds, Dotau, A, E = TILT(checkie, 'euclidean', init_points)

plt.imshow(A)
plt.show()
