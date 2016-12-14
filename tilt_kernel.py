from jacobi import *
from inner_IALM_constraints import *
from utils import *


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
    focus_size = np.floor(focus_size)

    image_size = input_image.shape
    focus_center = np.zeros((2,1))
    focus_center[0] = int((focus_size[1])/2)
    focus_center[1] =int((focus_size[0])/2)
    A_scale = 1

    UData = [1-image_center[0]-1, image_size[1]-image_center[0]-1]
    VData = [1-image_center[1]-1, image_size[0]-image_center[1]-1]
    XData = [1-focus_center[0]-1, focus_size[1]-focus_center[0]-1]
    YData = [1-focus_center[1]-1, focus_size[0]-focus_center[1]-1]

    #inp_im = np.hstack((np.zeros((50, 1)),input_image, np.zeros((50,1))))
    #inp_im = np.vstack((np.zeros((1, 52)),inp_im, np.zeros((1, 52))))
    #input_image = input_image.astype(np.uint8)

    #input_du = (inp_im[2:,:] - inp_im[:-2,:])[:,1:-1]
    #input_dv = (inp_im[:,2:] - inp_im[:,:-2])[1:-1, :]
    input_du = sobel(input_image, 1)
    input_dv = sobel(input_image, 0)

    #sobel
    #input_du
    Dotau_series = []

    tfm_matrix=initial_tfm_matrix
    #tfm=fliptform(maketform('projective', tfm_matrix'));
    ##CREATE INVERSE
    tfm_matrix = initial_tfm_matrix.T
    trfrm = transform.SimilarityTransform(tfm_matrix)
    #Dotau = transform.warp(input_image, trfrm)
    Dotau = transform_rak(input_image, tfm_matrix)

    #print "dotau: {0}".format(np.sum(Dotau))

    Dotau_series.append(Dotau)
    initial_image = Dotau

    #du = transform.warp(input_du, trfrm)
    #dv = transform.warp(input_dv, trfrm)
    du = transform_rak(input_du, tfm_matrix)
    dv = transform_rak(input_dv, tfm_matrix)


    #plt.imshow(du)
    #plt.show()

    du= du / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*du))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau
    dv= dv / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*dv))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau

    #plt.imshow(du)
    #plt.show()
    A_scale = np.linalg.norm(Dotau, 'fro')
    Dotau = Dotau.astype(float) / np.linalg.norm(Dotau, 'fro')
    #print "dotau_upd: {0}".format(np.sum(Dotau))
    tau = tfm2para(tfm_matrix, XData, YData, mode)
    J = jacobi(du, dv, XData, YData, tau, mode)
    S = constraints(tau, XData, YData, mode)

    print J[0, :]
    print S
    outer_round=0
    pre_f=0

    while 1:
        outer_round=outer_round+1
        A, E, delta_tau = inner_IALM_constraints(Dotau, J, S)
        #plt.imshow(A)
        #plt.show()
        #print('asas')
        #if error_sign == -1:
        #    return
        #end

        #if mod(outer_round, outer_display_period) == 0:
        #    disp(['outer_round ',num2str(outer_round),  ', rank(A)=', num2str(rank(A)), ', ||E||_1=', num2str(sum(sum(abs(E))))]);
        #print 'tau dtau', tau, delta_tau

        tau=tau + delta_tau
        print "tau", tau
        print 'dtau', delta_tau
        print outer_round

        tfm_matrix=para2tfm(tau, XData, YData, mode).T
        #tfm=fliptform(maketform('projective', tfm_matrix'))
        trfrm = transform.SimilarityTransform(tfm_matrix)
        Dotau = transform_rak(input_image, tfm_matrix)

        #print "tfm mat\n", tfm_matrix
        #print "dotau_upd: {0}".format(np.sum(Dotau))

        Dotau_series.append(Dotau)
        # judge convergence
        if (outer_round >= outer_max_iter): #or (abs(f-pre_f) < outer_tol):UPDATE STOPCOND
            break

        #pre_f=f;

        du = transform_rak(input_du, tfm_matrix)
        dv = transform_rak(input_dv, tfm_matrix)

        du= du / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*du))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau
        dv= dv / np.linalg.norm(Dotau, 'fro') - (sum(sum(Dotau*dv))) / (np.linalg.norm(Dotau, 'fro')**3) * Dotau
        A_scale = np.linalg.norm(Dotau, 'fro')
        Dotau = Dotau / np.linalg.norm(Dotau, 'fro')

        #tau = tfm2para(tfm_matrix, XData, YData, mode)
        J = jacobi(du, dv, XData, YData, tau, mode)
        S = constraints(tau, XData, YData, mode)



    #Dotau=imtransform(input_image, tfm, 'bilinear', 'UData', UData, 'VData',
    #VData, 'XData', XData, 'YData', YData, 'size', focus_size);
    #Dotau_series{1}=Dotau;
    #initial_image=Dotau;
    return Dotau, A, E, tfm_matrix, UData, VData, XData, YData, A_scale, Dotau_series