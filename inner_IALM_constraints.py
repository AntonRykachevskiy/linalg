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

from utils import *

def inner_IALM_constraints(Dotau, J, S_J, tol=1e-7, c=1., mu=-5., max_iter=100):
    #well shit then I guess I need to check stuff mu=1.25/np.linalg.norm(Dotau)
    """
    inner_IALM_constraints will solve the programming:
    min ||A||_*+lambda*||E||_1   s.t.     Dotau+J*delta_tau=A+E and
                                          S*delta_tau=0
    via the Inexact ALM method.
    ---------------------------------input----------------------------------
    Dotau:            m-by-n image matrix.
    J:                m-by-n-by-p tensor.
    S_J:              c-by-p matrix.
    inner_para:       parameters.
    --------------------------------output----------------------------------
    A:                m-by-n matrix, low-rank part.
    E:                m-by-n matrix, error part.
    delta_tau:        step of tau.
    f:                objective funtion value.
    """

    #prep data
    #sorry for a bit of bydlokod
    #if mu == -5.:
    mu=1.25/np.linalg.norm(Dotau, 2)

    m, n = Dotau.shape
    E = np.zeros((m, n))
    A = np.zeros((m, n))
    p = J.shape[2]
    delta_tau = np.zeros((p,1))

    J_vec = np.reshape(J, (m*n,1))#check whether reshape is the same
    Jo = J_vec
    #print 'J-vec before', J_vec.shape, S_J.shape,

    #J_vec = np.asarray([J_vec, S_J]) #plz check here as well, prob use smth else IS THIS TRUE I AM NOT SURE
    J_vec = np.vstack((J_vec, S_J))

    #print J_vec[:10]

    #print J
    #print 'Together',
    pinv_J_vec = np.linalg.pinv(J_vec)
    inner_round = 0
    rho = 1.25
    lmbda = c / np.sqrt(m)

    Y_1 = Dotau
    #print "Y_1:\n:", Y_1
    norm_two = np.linalg.norm(Y_1, 2)
    norm_inf = np.linalg.norm(Y_1.reshape(m*n, 1), np.inf)/lmbda #there was some mystery, what does norm_inf=norm(Y_1(:), inf)/lambda mean
    dual_norm = max(norm_two, norm_inf)
    Y_1 = Y_1 / dual_norm;
    Y_2 = np.zeros((S_J.shape[0], 1)) #this needs to be checked
    d_norm = np.linalg.norm(Dotau, 'fro')
    error_sign = 0
    #first_f = np.sum(np.linalg.svd(Dotau)) #this is so weird and I don't understand its purpose
    #stuff stuff stuff stuff stuff


    #print "2", norm_two
    #print "inf", norm_inf
    #print "dual", dual_norm
    #print "d", d_norm

    #let's now see what's in the main loop
    #ill just initialize stop_criterion to 10tol or something
    stop_criterion = 10*tol
    #max_iter = 50
    while (stop_criterion > tol) and (inner_round < max_iter):

        inner_round += 1
        temp_0 = Dotau + np.reshape(np.dot(Jo, delta_tau), (m,n)) + Y_1 / mu

        temp_1 = temp_0 - E
        U, S, V = np.linalg.svd(temp_1, full_matrices=False)
        S = np.diag(S)
        #print "U", U
        #print "S", S
        #print "V", V
        shrinkage_S =(S > 1/mu).astype(int) * (S - 1/mu) #because S is non-negative
        A = U.dot(shrinkage_S).dot(V)


        #print A.shape
        temp_2 = temp_0 - A
        E = (temp_2 > lmbda/mu).astype(int) * (temp_2 - lmbda/mu) + (temp_2 < -lmbda/mu).astype(int) * (temp_2 + lmbda/mu)
        #here is some MYSTERIOUS f, I do not understand its purpose
        temp_3 = A + E - Dotau - Y_1 / mu
        temp_3 = np.reshape(temp_3, (m*n,1))
        #temp_3 = np.asarray([temp_3, -Y_2/mu]) #CHECK CHECK CHECK CHECK
        temp_3 = np.vstack((temp_3, -Y_2/mu))

        #print 'ialm',pinv_J_vec,temp_3
        delta_tau = pinv_J_vec.dot(temp_3)
        derivative_Y_1 = Dotau - A - E + np.reshape(Jo.dot(delta_tau), (m, n))
        derivative_Y_2 = S_J.dot(delta_tau)
        Y_1 = Y_1 + derivative_Y_1 * mu  #mu is scalar
        Y_2 = Y_2 + derivative_Y_2 * mu

        stop_criterion=np.sqrt(np.linalg.norm(derivative_Y_1, 'fro')**2 + np.linalg.norm(derivative_Y_2, 2)**2)/d_norm
        mu=mu*rho

        #there is some mysterious error-catching, I didn't understand that part

    return A, E, delta_tau