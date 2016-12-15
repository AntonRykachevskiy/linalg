#prepare data for lowest resolution, whatever it means
import numpy as np
total_scale=np.floor(np.log2(np.min(args['focus_size']/args['focus_threshold']))) #whatever it means
downsample_matrix=np.array([[0.5, 0, 0], [0, 0.5, 0], [0, 0, 1]])
scale_matrix = np.linalg.matrix_power(downsample_matrix,total_scale)
tfm=maketform('projective', np.transpose(scale_matrix))############## don't know the hell is this
if args['blur'] == 1:
    input_image = GaussianBlur(args['input_image'], (5,5), 4)
else:
    input_image = args['input_image']
input_image=imtransform(input_image, tfm, 'bicubic') ##### you know better how to do it
if np.shape(input_image)[2]:
    input_image = rgb2gray(input_image) ##### you know how to handle it
input_image = input_image.astype(float)
initial_tfm_matrix = np.dot(np.dot(scale_matrix , args.initial_tfm_matrix),np.inv(scale_matrix))
center = floor(transform_point(args.center, scale_matrix))#### mad crazy functions
focus_size = np.floor(args['focus_size'] / (2**total_scale))
f_branch = zeros(3, 2 * args['branch_accuracy'] + 1)
#Dotau_branch = cell(3, 2 * args['branch_accuracy'] + 1)### array of matrices, don't exact equivalent in
Dotau_branch = np.empty((3,2 * args['branch_accuracy'] + 1,3,3))
result_tfm_matrix = np.empty((3,2 * args['branch_accuracy'] + 1,3,3))
    #step 2: design branch-and-bound method.
if args['mode'] == 'euclidian':
    max_rotation = args['branch_max_rotation']
    level = 1
    candidate_matrix = np.empty(1, 2 * args.branch_accuracy + 1,3,3)
    for i in range(2 * args.branch_accuracy + 1):
        candidate_matrix[0,i,:,:] = np.eye(3)
        theta = -max_rotation + (i - 1) * max_rotation / args['branch_accuracy']
        candidate_matrix[0,i,:2,:2]=np.array([[np.cos(theta), - np.sin(theta)],[np.sin(theta), np.cos(theta)]])
if args['mode'] == 'affine' or args['mode']:
    max_rotation = args['branch_max_rotation']
    max_skew = args['branch_max_skew']
    level = 3
    candidate_matrix = np.empty(3, 2 * args.branch_accuracy + 1, 3, 3)
    for i in range(2 * args.branch_accuracy + 1):
        candidate_matrix[0, i, :, :] = np.eye(3)
        theta = -max_rotation + (i - 1) * max_rotation / args['branch_accuracy']
        candidate_matrix[0, i, :2, :2] = np.array([[np.cos(theta), - np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        candidate_matrix[1,i,:,:] = np.eye(3)
        candidate_matrix[1, i,0, 1] = -max_skew + (i - 1) * max_skew / args['branch_accuracy']
        candidate_matrix[2,i,:,:] = np.eye(3)
        candidate_matrix[2, i, 1, 0] = -max_skew + (i - 1) * max_skew / args['branch_accuracy']
gap = 5
BLACK_MATRIX = np.zeros(focus_size[0] * level + gap * (level - 1), focus_size[1] * (2 * args.branch_accuracy + 1) + gap * 2 * args['branch_accuracy'])
    #step 3  - begin branch and bound
normal_outer_max_iter = args['outer_max_iter']
normal_display_inter = args['display_result']
args['outer_max_iter'] = 1#for debug, set it to 1;
args['display_result'] = 0
for i in range(level):
    for j in range (2 * args.branch_accuracy + 1):
        tfm_matrix = np.inv(candidate_matrix[i, j,:,:] * inv(initial_tfm_matrix))### don't want mess with this
        args['figure_no'] = (i - 1) * level + j
        args['save_path'] = [] ### some crazy strange stuff
        image_size = np.shape(args['input_image'])
        image_center = np.floor(args['center'])
        focus_center = np.zeros((2, 1))
        focus_center[0] = np.floor((1 + args['focus_size'][1]) / 2)
        focus_center[1] = np.floor((1 + args['focus_size'][0]) / 2)
        UData = np.array([1 - image_center[0], image_size[1] - image_center[0]])
        VData = np.array([1 - image_center[1], image_size[0] - image_center[1]])
        XData = np.array([1 - focus_center[0], args['focus_size'][1] - focus_center[0]])
        YData = np.array([1 - focus_center[1], args['focus_size'][0] - focus_center[1]])
        tfm = fliptform(maketform('projective', np.transpose(tfm_matrix))) ### I don't know how to use this function
        Dotau = imtransform(input_image, tfm, 'bilinear', 'XData', XData, 'YData', YData, 'UData', UData, 'VData',VData, 'Size', focus_size) ###the same
        Dotau = Dotau / np.linalg.norm(Dotau, 'fro')
        U, S, V = np.linalg.svd(Dotau)
        f = np.sum(S) ##### should be carefull
        start = [(focus_size[0] + gap) * (i - 1) + 1, (focus_size[1] + gap) * (j - 1) + 1]
        BLACK_MATRIX[start[0]:(start[0] + focus_size[1] - 1), start[1]:(start[1] + focus_size[1] - 1)] = Dotau ## don't know what it is
        f_branch[i, j] = f
        Dotau_branch[i, j,:,:] = Dotau
        result_tfm_matrix[i, j,:,:] = tfm_matrix
    index = np.argmin(f_branch[i,:])
    value = np.amin(f_branch[i,:])
    initial_tfm_matrix = result_tfm_matrix[i, index,:,:] ## I don't know where i goes
#step 4
initial_tfm_matrix = np.dot(np.dot(np.inv(scale_matrix),initial_tfm_matrix),scale_matrix)
args['outer_max_iter']=normal_outer_max_iter
args['display_result']=normal_display_inter


