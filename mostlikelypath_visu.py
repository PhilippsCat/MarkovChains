import numpy as np
def mostlikelypath_visu(start, T, pathlength):
    M = np.zeros((pathlength,T.shape[0]))
    M[0,:] = start
    dispmatrix = M
    for i in range(1,pathlength):
        M[i,:] = dot(T.transpose(),start)
        dispmatrix[i,:] = M[i,:]*M[i,:]*M[i,:]*M[i,:]
        dispmatrix[i,:] = dispmatrix[i,:] / sum(dispmatrix[i,:])
        T = dot(T,T) 
    plt.imshow(dispmatrix.transpose())
    return
