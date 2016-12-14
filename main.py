from TILT import *

if __name__ == '__main__':
    check = io.imread('cm_1.jpg')
    print np.amax(check)

    init_points = np.asarray([[0, check.shape[1]], [0, check.shape[0]]])

    print "shape",  check.shape

    checkie = np.zeros((check.shape[0],check.shape[1],3))
    checkie[:,:,0] = check
    checkie[:,:,1] = check
    checkie[:,:,2] = check
    Ds, Dotau, A, E = TILT(checkie, 'euclidean', init_points)


    plt.imshow(A)
    #plt.imshow(E)
    plt.show()


    plt.imshow(E)
    #plt.imshow(E)
    plt.show()
