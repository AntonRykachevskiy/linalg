from jacobi import *
from inner_IALM_constraints import *
from utils import *
import cv2


def polina_transform(input_image, tfm_matrix, UData, VData, XData, YData):

    final = np.eye(3)

    #first, translate to a new center:
    u_trans = np.eye(3)
    u_trans[0,2] = UData[0]
    u_trans[1,2] = VData[0]

    #then apply the transform
    M = np.eye(3)
    M = np.linalg.inv(tfm_matrix)

    #then do a transform according to XData and YData
    x_trans = np.eye(3)
    x_trans[0,2] = -XData[0]
    x_trans[1,2] = -YData[0]

    final = x_trans.dot(M).dot(u_trans)

    img_1 =  cv2.warpPerspective(input_image, final, (np.sum(np.abs(XData)).astype(int),np.sum(np.abs(YData)).astype(int)))

    return img_1


def transform_rak(image, tfm_matrix):
    shift_y, shift_x = np.array(image.shape[:2]) / 2.
    tf_rotate = transform.SimilarityTransform(tfm_matrix)
    tf_shift = transform.SimilarityTransform(translation=[-shift_x, -shift_y])
    tf_shift_inv = transform.SimilarityTransform(translation=[shift_x, shift_y])

    image_rotated = transform.warp(image, (tf_shift + (tf_rotate + tf_shift_inv)).inverse)

    return image_rotated


def tilt_kernel(input_image, mode, center, focus_size, initial_tfm_matrix, para):
    outer_tol = 5e-5
    outer_max_iter = 100
    outer_display_perioud = 1


    if input_image.shape[2] > 1:
        input_image = input_image[:,:,0]*0.299 + input_image[:,:,1]*0.587 + input_image[:,:,2]*0.144

    input_image = input_image.astype(float)

    image_center = np.floor(center)
    print 'im_c', image_center
    fointcus_size = np.floor(focus_size)
    print 'fs', focus_size
    image_size = input_image.shape

    print image_size
    focus_center = np.zeros((2,1))
    focus_center[0] = int((focus_size[1])/2)
    focus_center[1] =int((focus_size[0])/2)
    A_scale = 1

    UData = [1-image_center[0], image_size[1]-image_center[0]-1]
    VData = [1-image_center[1], image_size[0]-image_center[1]-1]
    XData = [1-focus_center[0], focus_size[1]-focus_center[0]-1]
    YData = [1-focus_center[1], focus_size[0]-focus_center[1]-1]

    #inp_im = np.hstack((np.zeros((input_image.shape[0], 1)),input_image, np.zeros((input_image.shape[0],1))))
    #inp_im = np.vstack((np.zeros((1, input_image.shape[1] +2 )),inp_im, np.zeros((1, input_image.shape[1]  + 2))))
    #input_image = input_image.astype(np.uint8)

    #input_du = (inp_im[2:,:] - inp_im[:-2,:])[:,1:-1]
    #input_dv = (inp_im[:,2:] - inp_im[:,:-2])[1:-1, :]
    input_du = sobel(input_image, 1)
    input_dv = sobel(input_image, 0)

    Dotau_series = []

    tfm_matrix=initial_tfm_matrix
    #Piece of rotation code
    #tfm_matrix = initial_tfm_matrix.T
    #Dotau = transform_rak(input_image, tfm_matrix)

    Dotau = polina_transform(input_image, tfm_matrix, UData, VData, XData, YData)

    #print "dotau: {0}".format(np.sum(Dotau))

    Dotau_series.append(Dotau)

    #du = transform_rak(input_du, tfm_matrix)
    #dv = transform_rak(input_dv, tfm_matrix)
    du = polina_transform(input_du, tfm_matrix, UData, VData, XData, YData)
    dv = polina_transform(input_dv, tfm_matrix, UData, VData, XData, YData)

    du= du / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*du))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau
    dv= dv / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*dv))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau

    plt.imshow(du, cmap='gray')
    plt.show()

    A_scale = np.linalg.norm(Dotau, 'fro')
    Dotau = Dotau.astype(float) / np.linalg.norm(Dotau, 'fro')
    print 'max', np.max(du)
    print np.max(dv)
    print np.min(dv)

    tau = tfm2para(tfm_matrix, XData, YData, mode)
    print tfm_matrix
    print XData
    print YData

    print 'ttttt', tau
    J = jacobi(du, dv, XData, YData, tau, mode)
    S = constraints(tau, XData, YData, mode)

    print S

    outer_round=0
    pre_f=0

    print tau.shape

    while 1:
        outer_round=outer_round+1
        A, E, delta_tau = inner_IALM_constraints(Dotau, J, S)

        tau=tau + delta_tau
        print "tau", tau
        print 'dtau', delta_tau
        print outer_round

        tfm_matrix=para2tfm(tau, XData, YData, mode)

        #Dotau = transform_rak(input_image, tfm_matrix)
        Dotau = polina_transform(input_image, tfm_matrix, UData, VData, XData, YData)


        Dotau_series.append(Dotau)
        # judge convergence
        if (outer_round >= outer_max_iter):
            break

        #du = transform_rak(input_du, tfm_matrix)
        #dv = transform_rak(input_dv, tfm_matrix)
        du = polina_transform(input_du, tfm_matrix, UData, VData, XData, YData)
        dv = polina_transform(input_dv, tfm_matrix, UData, VData, XData, YData)

        du= du / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*du))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau
        dv= dv / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*dv))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau
        A_scale = np.linalg.norm(Dotau, 'fro')
        Dotau = Dotau / np.linalg.norm(Dotau, 'fro')

        J = jacobi(du, dv, XData, YData, tau, mode)
        S = constraints(tau, XData, YData, mode)


    return Dotau, A, E, tfm_matrix, UData, VData, XData, YData, A_scale, Dotau_series