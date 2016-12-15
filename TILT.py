from cv2 import GaussianBlur

from tilt_kernel import *

def TILT(input_image, mode, init_points, **kwargs):
    args = parse_TILT_arguments(**kwargs)
    args['mode'] = mode
    args['input_image'] = input_image
    args['initial_points'] = np.floor(init_points)
    args['focus_size'] = np.asarray([init_points[1,1]-args['initial_points'][1,0], args['initial_points'][0,1]-args['initial_points'][0,0]])
    args['center'] = np.floor(np.mean(init_points, axis=1))
    args['focus_center'] = np.floor(np.asarray([args['focus_size'][1],  args['focus_size'][0]])/2.) #onwes were removed
    focus_center = args['focus_center']
    XData = np.asarray([1-focus_center[0], args['focus_size'][1]-focus_center[0]])
    YData = np.asarray([1-focus_center[1], args['focus_size'][0]-focus_center[1]])
    args['pyramid'] = 1
    #now, afterwards there is nasty stuff about compute homography
    #hopefully, I can avoid it and initialize tfm to identity

    original_args = args

    #creating boundaries of the image
    expand_rate = 0.8 #WTF is this shit
    initial_points = args['initial_points']
    left_bound = np.ceil(max(initial_points[0,0] - expand_rate * (initial_points[0,1] - initial_points[0,0]), 0))
    right_bound = np.floor(min(initial_points[0,1] + expand_rate * (initial_points[0,1] - initial_points[0,0]), input_image.shape[1]));
    top_bound = np.ceil(max(initial_points[1,0] - expand_rate * (initial_points[1,1] - initial_points[1,0]), 0))
    bottom_bound = np.floor(min(initial_points[1,1] + expand_rate*(initial_points[1,1] - initial_points[1,0]), input_image.shape[0]));
    new_image = np.zeros((bottom_bound - top_bound , right_bound - left_bound , input_image.shape[2]))

    print 'bounds'
    print left_bound
    print right_bound
    print  top_bound
    print bottom_bound

    for c in range(input_image.shape[2]):
        new_image[:,:,c]=input_image[top_bound:bottom_bound, left_bound:right_bound, c] #maybe miss one pixel? but whatever

    #new_image = input_image #TEMP

    #note - here I didn't copy this: args.input_image=uint8(new_image),
    #because that would be handled at input by scikit
    #sorry, to tired to make it not bydlocode
    #args['input_image'] = img_as_ubyte(new_image) #hope this is right
    args['input_image'] = new_image


    args['center'] = args['center'] + np.asarray([1-left_bound, 1-top_bound])


    pre_scale_matrix=np.eye(3)
    #here should be the part "down-sample the image if the focus is too large"
    #but guess what
    #it won't be too large
    #at least for initial testing

    min_length = original_args['focus_size']

    initial_tfm_matrix = args['initial_tfm_matrix']
    initial_tfm_matrix = np.linalg.inv(pre_scale_matrix).dot(initial_tfm_matrix).dot(pre_scale_matrix)
    args['initial_tfm_matrix'] = initial_tfm_matrix
    args['focus_size'] = np.around(args['focus_size']/pre_scale_matrix[0,0])
    args['center'] = args['center']/pre_scale_matrix[0,0]
    parent_path='./'





    if args['pyramid'] == 1:
        downsample_matrix = np.array([[0.5, 0, 0], [0, 0.5, 0], [0, 0, 1]])
        upsample_matrix = np.linalg.inv(downsample_matrix)
        total_scale = np.ceil(np.max(np.log2(np.min(args['focus_size'])/args['focus_threshold']), 0));

        for scale in range(int(total_scale),-1,-1):
            # begin each level of the pyramid
            if total_scale - scale >= args['pyramid_max_level']:
                break

            # Blur if required
            #BULA HERNIA
            if args['blur'] == 1 and scale == 0:
                input_image = GaussianBlur(args['input_image'], (5,5), 4) ##GIVE KERNEL
                #input_image=imfilter(args.input_image, fspecial('gaussian', ceil(args.blur_kernel_size_k*2^scale), ceil(args.blur_kernel_sigma_k*2^scale)));
            else:
                input_image = args['input_image']


            # prepare image and initial tfm_matrix
            scale_matrix = np.linalg.matrix_power(downsample_matrix, scale)
            #tfm = maketform('projective', scale_matrix');
            #input_image=imtransform(input_image, tfm, 'bicubic');

            input_image = polina_transform(input_image, scale_matrix, inv_flag = False)

            print scale_matrix
            tfm_matrix = scale_matrix.dot(initial_tfm_matrix).dot(np.linalg.inv(scale_matrix))

            center = np.floor(transform_point(args['center'], scale_matrix));

            focus_size=np.floor(args['focus_size']/(2 ** scale))

            #args['save_path'] = fullfile(parent_path, ['pyramid', num2str(scale)]);
            #args.figure_no=100+total_scale-scale+1;
            Dotau, A, E, tfm_matrix, UData, VData, XData, YData, A_scale, Dotau_series = tilt_kernel(input_image, args['mode'], center, focus_size, tfm_matrix, args)
            # update tfm_matrix of the highest-resolution level.
            initial_tfm_matrix = np.linalg.inv(scale_matrix).dot(tfm_matrix).dot(scale_matrix)
            args['outer_tol']=args['outer_tol']*args['outer_tol_step']

        tfm_matrix=initial_tfm_matrix

    else:
        img_size=np.shape(args['input_image'])
        img_size=img_size[:2]
        img_size=img_size[:2]
        input_image = GaussianBlur(args['input_image'], (5,5), 4)

        args['figure_no'] = 101
        args['save_path'] = os.path.join(parent_path, 'some_name')
        #print args['center']
        #print args['focus_size']
        Dotau, A, E, tfm_matrix, UData, VData, XData, YData, A_scale, Dotau_series = tilt_kernel(args['input_image'], args['mode'], args['center'], args['focus_size'], args['initial_tfm_matrix'], args)

    args = original_args
    tfm_matrix = np.dot(pre_scale_matrix,np.dot(tfm_matrix,np.linalg.inv(pre_scale_matrix)))

    focus_size=args['focus_size']
    image_size=np.shape(args['input_image'])
    image_size=image_size[:2]
    image_center=args['center']
    focus_center=np.zeros((2, 1))
    focus_center[0]=np.floor((1+args['focus_size'][1])/2)
    focus_center[1]=np.floor((1+args['focus_size'][0])/2)

    UData=np.array([1-image_center[0], image_size[1]-image_center[0]])
    VData=np.array([1-image_center[1], image_size[0]-image_center[1]])
    XData=np.array([1-focus_center[0], args['focus_size'][1]-focus_center[0]])
    YData=np.array([1-focus_center[1], args['focus_size'][0]-focus_center[1]])

    #display the result    #I NEED TO MAKE THE SAME CORRECTION BELOW - P.
    #if args['display_result'] == 1:
       # plt.figure(99)
        #plt.imshow(args[input_image], [], 'DisplayRange', [0, max(max(max(args.input_image)))])
        # I don't what are parameters are here, may be we should omit them
        #X1=[args['initial_points'][0,0], args['initial_points'][0,0], args['initial_points'][0,1], args['initial_points'][0,1], args['initial_points'][0,0]]
        #Y1=[args['initial_points'][1,0], args['initial_points'][1,1], args['initial_points'][1,1], args['initial_points'][1,0], args['initial_points'][1,0]]
        #plt.hold(True)
        #plt.plot(X1, Y1, 'r-')
        #pt_top_left = transform_point([XData[0], YData[0]], tfm_matrix)+image_center
        #pt_bottom_left = transform_point([XData[0], YData[1]], tfm_matrix)+image_center
        #pt_bottom_right = transform_point([XData[1], YData[1]], tfm_matrix)+image_center
        #pt_top_right = transform_point([XData[1], YData(1)], tfm_matrix)+image_center
        #X2=[pt_top_left[0], pt_bottom_left[0], pt_bottom_right[0], pt_top_right[0], pt_top_left[0]]
        #Y2=[pt_top_left[1], pt_bottom_left[1], pt_bottom_right[1], pt_top_right[1], pt_top_left[1]]
        #plt.plot(X2, Y2, 'g-')

    return  Dotau_series, Dotau, A, E