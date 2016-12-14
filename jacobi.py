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


def jacobi(du, dv, XData, YData, tau, mode):
    """ jacobi() will calculate the jacobi matrix.
    ------------------------------input--------------------------------------
    du:           m-by-n matrix, derivative along x-axis;
    dv:           m-by-n matrtix, derivative along y-axis;
    UData:        1-by-2 vector, X-range of the image.
    VData:        1-by-2 vector, Y-range of the image.
    tfm_matrix:   3-by-3 matrix, transformation matrix.
    mode:         one of 'affine', 'affine_notranslation', 'homography',
    'homography_notranslation'
    ------------------------------output-------------------------------------
    J:            m-by-n-by-p tensor, jacobi matrix.

    """
    [m, n] = du.shape
    [X0, Y0] = np.meshgrid(range(XData[0],XData[1]+1), range(YData[0],YData[1]+1))

    if mode == 'euclidean':
        tau = [tau]
        J = np.zeros((m, n, 1))
        J[:,:,0] = du * (- X0 * np.sin(tau[0]) - Y0 * np.cos(tau[0])) + dv * (X0 * np.cos(tau[0]) - Y0 * np.sin(tau[0]))
    #if mode == 'affine':

    return J