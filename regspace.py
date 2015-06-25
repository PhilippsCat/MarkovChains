def regspace(X , dmin):

    T = X.shape[1]
    Y = list()
    Y.append(X[:,0])
    
    for t in range(2,T):
        if ( np.min ( np.linalg.norm(X[:,t] - Y)) < dmin):
            Y.append(X[:,t])
    return np.asarray(Y)