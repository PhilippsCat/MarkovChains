import numpy as np

def regspace(X , dmin):


    
    # Reshape the data 
    if (len(X.shape) == 2):
        len_X = max( X.shape[0], X.shape[1])
    else:
        len_X = len(X)
    X = X.reshape(-1,len_X)
    T = len_X
    Y = list()
    Y.append(X[:,0])
    
    # Do the regular space algorithm
    for t in range(1,T):
        if ( np.min ( np.absolute(X[:,t]*np.ones(len(Y)) - Y)) > dmin):
            Y.append(X[:,t])
     
    # The format of Y is changed from list to np.array
    Y = np.asarray(Y)

    
    closest_centers = np.zeros(T)
    # values are mapped to their closest centers
    for i in range(0,T):
        closest_center_index = np.argmin( np.absolute(X[:,i]*np.ones(len(Y)) - Y.transpose()))
        closest_centers[i] = Y[closest_center_index]
    
    print closest_centers
    # 'uniques' is the list of states that occur in 'closest_centers'
    uniques = list(set(closest_centers))
    # we map the values in 'closest_centers' to integer values
    int_states = np.zeros(T)
    for j  in range(0,T):
        int_states[j] = uniques.index(closest_centers[j])

    
    
    return int_states, Y
