import numpy as np

def kmeans(data, k):
    """
    Performs a k-means clustering on a given data set, which should have the
    dimension d x n, where d is the dimension of each data point and n is the
    total number of data points. k is the number of clusters.
    Usage: kmeans(data, k)
    Return: Array that maps each data point to a cluster, array which
    contains the centers of the clusters.
    """
    n = np.shape(data)[1] #suppose data steps are in y direction  
    means = data[:,np.random.choice(n,k)]

    def closest_means(current_means):
        associated_mean = np.ones(n)
        for i in range(0,n):
            dists = np.linalg.norm(np.vstack(data[:,i])-current_means, axis=0)
            associated_mean[i] = np.argmin(dists)
        return associated_mean
    
    def calculate_centroids(means):
        for i in range(0,k):
                # if there is no point associated to the current mean,
                #don't change anything
                if np.sum(point_to_means==i)==0:        
                    continue
                means[:,i]=np.sum(data[:, point_to_means==i],1)/np.sum(point_to_means==i)
        return means

    changed = True
    ptm_old = np.zeros(n)
    while changed: 
        point_to_means = closest_means(means)            
        if np.array_equal(ptm_old, point_to_means):
            changed = False
        ptm_old = point_to_means
        means = calculate_centroids(means)
    
    return point_to_means, means
