import matplotlib.pyplot as plt
from skimage import data, io, filters



from TILT import *

if __name__ == '__main__':
    check = cv2.imread('building_.jpg')


    #init_points = np.asarray([[0, check.shape[1]], [0, check.shape[0]]])

    init_points = np.asarray([[40, 80], [40, 80]])

    plt.imshow(check[init_points[0][0]: init_points[0][1], init_points[1][0]: init_points[1][1]])
    plt.show()


    #checkie = np.zeros((check.shape[0],check.shape[1],3))
    #checkie[:,:,0] = check
    #checkie[:,:,1] = check
    #checkie[:,:,2] = check
    Ds, Dotau, A, E = TILT(check, 'homography', init_points)


    plt.imshow(A)
    #plt.imshow(E)
    plt.show()


    plt.imshow(E)
    #plt.imshow(E)
    plt.show()

    print np.sum(E)
