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
