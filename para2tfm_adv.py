import numpy as np
def para2tfm(tau, XData, YData, mode):
    tfm_matrix = np.identity(3)
    if mode == 'euclidian':
        tfm_matrix[:2,:2]=np.array(([np.cos(tau[0]),-np.sin(tau[0])],[np.sin(tau[0]),np.cos(tau[0])]))
    if mode == 'affine':
        tfm_matrix[0,:2] = np.transpose(tau[:2]) #should accurate here
        tfm_matrix[1,:2] = np.transpose(tau[2:])
    if mode == 'homography':
        X = [XData[0], XData[1], XData[1], XData[0]]
        Y = [YData[0], YData[0], YData[1], YData[1]]
        temp = np.reshpe(tau,(2,4))
        U = temp[0,:]
        V = temp[1,:]
        A = np.zeros((8, 8))
        b = np.zeros((8, 1))
        insert_A = np.zeros((2,8))
        insert_b = np.zeros((2,1))
        for i in range(4):
            insert_A[0,:]=[0, 0, 0, - X[i], - Y[i], - 1, V[i] * X[i], V[i] * Y[i]]
            insert_b[0] = -V[i]
            insert_A[1,:]=[X[i], Y[i], 1, 0, 0, 0, - U[i] * X[i], - U[i] * Y[i]]
            insert_b[i] = U[i]
            A[2 * i +1:2 * i,:] = insert_A #error might occur here
            b[2 * i + 1:2 * i] = insert_b
        solution = A/b #don't know what happening here
    tfm_matrix = np.transpose(np.reshape([solution,1],(3,3)))




















