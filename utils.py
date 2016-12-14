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



def parse_TILT_arguments(**kwargs):
    my_args = kwargs

    if 'initial_tfm_matrix' not in my_args:
        my_args['initial_tfm_matrix'] = np.eye(3)

    if 'outer_tol' not in my_args:
        my_args['outer_tol'] = 1e-4

    if 'outer_max_iter' not in my_args:
        my_args['outer_max_iter'] = 50

    if 'outer_display_period' not in my_args:
        my_args['outer_display_period'] = 1

    if 'inner_tol' not in my_args:
        my_args['inner_tol'] = 1e-4

    if 'inner_c' not in my_args:
        my_args['inner_c'] = 1

    if 'inner_mu' not in my_args:
        my_args['inner_mu'] = []  #questionable

    if 'inner_display_period' not in my_args:
        my_args['inner_display_period'] = 100

    if 'inner_max_iter' not in my_args:
        my_args['inner_max_iter'] = np.inf

    if 'blur' not in my_args:
        my_args['blur'] = 1

    if 'branch' not in my_args:
        my_args['branch'] = 1

    if 'pyramid' not in my_args:
        my_args['pyramid'] = 1

    if 'focus_threshold' not in my_args:  #when doing pyramid the smallest focus_edge we can tolerate.
        my_args['focus_threshold'] = 50

    if 'outer_tol_step' not in my_args:  #as resolution goes high how relaxed should the outer_tol be
        my_args['outer_tol_step'] = 10

    if 'blur_kernel_size_k' not in my_args:  # neighbourhood scalar for the size of the blur kernel.
        my_args['blur_kernel_size_k'] = 2

    if 'blur_kernel_sigma_k' not in my_args: # standard derivation scalar for blur kernel.
        my_args['blur_kernel_sigma_k'] = 2

    if 'pyramid_max_level' not in my_args: # number of pyramid levels we want to act on.
        my_args['pyramid_max_level'] = 2

    if 'branch_max_iter' not in my_args: # in each branch, how much iteration we take.
        my_args['branch_max_iter'] = 10

    if 'branch_max_accuracy' not in my_args: # higher means smaller step-width.
        my_args['branch_max_accuracy'] = 5

    if 'display_result' not in my_args:
        my_args['display_result'] = 1

    if 'focus_size' not in my_args:
        my_args['focus_size'] = []

    if 'save_path' not in my_args:
        my_args['save_path'] = []

    #MORE STUFF FOR BRANCH AND BOUND SKIPPED HERE! ADD LATER
    #Stuff skipped here because I do it outside

    return my_args



def constraints(tau, XData, YData, mode):
    """constraints() will get the linearize constraints of tau according to mode.
    -----------------------------input--------------------------------------
    tau:          p-by-1 real vector.
    mode:         one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 'homography',
    'homography_notranslation'.
    ----------------------------output--------------------------------------
    linearized constraints on tau.

    """

    if mode == 'euclidean':
        S = np.zeros((2,1))

    return S
    #if mode == 'affine':



def para2tfm(tau, XData, YData, mode):
    """para2tfm will turn tau to tfm_matrix according to mode.
       ----------------------------input---------------------------------------
       tau:      p-by-1 vector
       mode:     one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 'homography',
                'homography_notranslation'
       ----------------------------output--------------------------------------
       tfm_matrix:   3-by-3 transform matrix.
    """
    tfm_matrix = np.eye(3)
    if mode == 'euclidean':
        tfm_matrix[0:2,0:2] = [[np.cos(tau[0]), -np.sin(tau[0])],
                               [np.sin(tau[0]), np.cos(tau[0])]]

    elif mode == 'affine':
        tfm_matrix[0,0:2] = tau[0:2].H
        tfm_matrix[1,0:2] = tau[2:].H

    else:
        print 'no param'

    return tfm_matrix


def tfm2para(tfm_matrix, XData, YData, mode):
    """tfm2para will transpose tfm_matrix to its corresponding parameter.
     -------------------------input------------------------------------------
     tfm_matrix:       3-by-3 matrix.
     mode:             one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation',
                       'homography', 'homography_notranslation'
     -------------------------output-----------------------------------------
     tau:              p-by-1 real vector.
    """
    if mode == 'euclidean':
        tau = np.arccos(tfm_matrix[0, 0])
        if tfm_matrix[1, 0] < 0:
            tau *= -1.
        tau = np.array(tau)

    elif mode =='affine':
        tau = np.zeros((4, 1))
        tau[0:2] = tfm_matrix[0,0:2].H
        tau[2:] = tfm_matrix[1,0:2].H

    else:
        print 'wut'

    return tau


def transform_point(input_pt, tfm_matrix):
    #if size(input_pt, 1)==1
    if input_pt.shape[0] == 1: #we make input into a column vector
        input_pt = input_pt.H
        b_row = 1
    else:
        b_row = 0

    input_pt = np.asarray([input_pt, 1]) #please match shapes please
    output_pt = tfm_matrix.dot(input_pt)
    output_pt = output_pt/output_pt[2]
    output_pt = output_pt[0:2, 0]
    if b_row == 1:
        output_pt = output.H

    return output_pt