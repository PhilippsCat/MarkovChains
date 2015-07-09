import numpy as np

def calculate_centroids(data, means, point_to_means):
    """
    Calculate the centers of sets of vectors of arbitrary length (all same).
    Input: 
        data: all vectors
        means: position of current centers
        point_to_means: maps vectors to associated centers
                        (data[i] corresponds to means[point_to_means[i]])
    Output: cluster centers of sets of vectors given by point_to_means
    """
    k = np.shape(means)[0]
    new_means = means
    for i in range(0,k):
            # if there is no point associated to the current mean,
            # don't change anything
            if np.sum(point_to_means==i)==0:        
                continue
            new_means[i] = np.sum(data[point_to_means==i],0)/np.sum(point_to_means==i)
    return new_means

def closest_means(data, means):
    """
    Creates the mapping point_to_means described in calculate_centroids
    """
    n = np.shape(data)[0] 
    k = np.shape(means)[0]
    assoc_means = np.ones(n)
    for i in range(0,n):
        dists = np.linalg.norm(np.tile(data[i],(k,1))-means,axis=1)
        assoc_means[i] = np.argmin(dists)
    return assoc_means
 
def kmeans(data, k):
    """
    Performs a k-means clustering on a given data set, which should have the
    dimension n x d, where d is the dimension of each data point and n is the
    total number of data points. k is the number of clusters.
    Return: Array that maps each data point to a cluster, array which
    contains the centers of the clusters.
    """
    n = np.shape(data)[0] 
    means = data[np.random.choice(n,k, replace=False)]
   
    ptm_old = np.zeros(n)
    while True: 
        point_to_means = closest_means(data, means)            
        if np.array_equal(ptm_old, point_to_means):
            break
        ptm_old = point_to_means
        means = calculate_centroids(data, means, point_to_means)
    
    return point_to_means, means
    
def regspace(X , dmin):


    #This is the regular space clustering alogirthm. It needs a one- or multi-dimensional numpy-array X and a scalar value dmin
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
        int_states[j] = (np.where(np.all(uniques==closest_centers[j,:],axis=1)))[0]


    
    
    return int_states, Y
