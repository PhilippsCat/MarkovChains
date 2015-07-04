import numpy as np
import matplotlib.pyplot as plt
def mostlikelypath_visu(start, T, pathlength):
    M = np.zeros((pathlength,T.shape[0]))
    M[0,:] = np.zeros((1,T.shape[0]))
    M[0,start] = 1
    dispmatrix = M
    for i in range(1,pathlength):
        M[i,:] = np.dot(T.transpose(),M[0,:])
        dispmatrix[i,:] = M[i,:]*M[i,:]*M[i,:]*M[i,:]
        dispmatrix[i,:] = dispmatrix[i,:] / sum(dispmatrix[i,:])
        T = np.dot(T,T) 
    plt.imshow(dispmatrix.transpose())
    plt.show()
    return
    
