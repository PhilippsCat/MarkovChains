import numpy as np

def kmeans(data, k):
    n = np.shape(data)[1] #suppose data steps are in y direction  
    means = data[:,np.random.choice(n,k)]

    def closest_means(current_means):
        associated_mean = np.ones(n)
        for i in range(0,n):
            dists = np.linalg.norm(np.vstack(data[:,i])-current_means, axis=0)
            associated_mean[i] = np.argmin(dists)
        return associated_mean
    
    def calculate_centroids():
        for i in range(0,k):
                means[:,i]=np.sum(data[:, point_to_means==i],1)/np.sum(point_to_means==i)

    changed = True
    ptm_old = np.zeros(n)
    j = 0
    while changed: 
        j=j+1
        if j>1000:
            print "1000 reached"
            break
        point_to_means = closest_means(means)            
        if np.array_equal(ptm_old, point_to_means):
            changed = False
        ptm_old = point_to_means
        means = calculate_centroids()
    
    return point_to_mean, means

