import numpy as np

def regspace(X , dmin):


    
    # Reshape the data. It is taken care of the three types of one-dimensional data: [x,y,z], [[x,y,z]] and [ [x], [y], [z]]
    if (len(X.shape) == 2):
        len_X = X.shape[0]
        if (X.shape[0] == 1 and X.shape[1] > 1):
            X = X[0]
            len_X = len(X)
            X = X.reshape(-1,len_X)
            X = X.reshape(-1,len(X))
    else:
        len_X = len(X)
        X = X.reshape(-1,len_X)
        X = X.reshape(-1,len(X))

    T = len_X
    Y = list()
    X = X.transpose()
    Y.append(X[:,0])
    # Do the regular space algorithm
    for t in range(1,T):
        if ( np.min ( np.linalg.norm(X[:,t]*np.ones((len(Y),Y[0].shape[0])) - Y, axis = 1)) > dmin):
            Y.append(X[:,t])
     
    # The format of Y is changed from list to np.array
    Y = np.asarray(Y)

    closest_centers = np.zeros((T,Y[0].shape[0]))
    
    # values are mapped to their closest centers
    for i in range(0,T):
        closest_center_index = np.argmin( np.linalg.norm(X[:,i]*np.ones((len(Y),Y[0].shape[0])) - Y, axis = 1))
        closest_centers[i,:] = Y[closest_center_index]
    
    # 'uniques' is the list of states that occur in 'closest_centers'   
    ncols = closest_centers.shape[1]
    dtype = closest_centers.dtype.descr * ncols
    struct = closest_centers.view(dtype)
    uniques = np.unique(struct)
    uniques = uniques.view(closest_centers.dtype).reshape(-1, ncols)
    
    # we map the values in 'closest_centers' to integer values
    int_states = np.zeros(T)
    for j  in range(0,T):
        na = closest_centers[j,:]
        int_states[j] = (np.where(np.all(uniques==closest_centers[j,:],axis=1)))[0]


    
    
    return int_states, Y
