import matplotlib.pyplot as plt
import numpy as np



def plot_chain_frequencies(chain):
    #This method sorts the states by their number of occurrences in 
    # a chain und plots the numbers of occurrences in descending order. It requires
    # a 1D-numpy-array
    states = list(set(chain))
    numofstates = len(states)
    count = np.zeros([numofstates,1])
    counter = 0
    
    for i in range(0,numofstates):
        count[counter] = (chain == states[i]).sum()
        counter = counter + 1
  
    plt.plot(sorted(count, reverse=True))
    plt.xlabel("Rank of state (by number of occurrences)")
    plt.ylabel("Number of occurrences")
    plt.axis([0, max(states), 0, max(count)])
    plt.show()
    return count
    
    

def rank_by_std(T):
    n = T.shape[0]
    variances = np.zeros(n)
    for i in range(0,n):
        variances[i] = np.std(T[i,:])
    return np.lexsort((range(0,n),variances)), variances

def plot_2d_with_clusters(data, clusters):
    """ Plots a scatter plot of 2d data where the points are colored
    as their clusters suggest.
    """
    plt.scatter(data[:,0],data[:,1], c=clusters[0])
    plt.show()
    
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
