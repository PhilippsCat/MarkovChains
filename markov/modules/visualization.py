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
# This method ranks the states by the variances that their respective rows
# in the transition matrices have. The last state in the output is the state with
# the highest variances between transition probabilities
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
# This method needs a given start distribution and plots the most probable paths
# that the induced process is going to take in the brightest colors.
    M = np.zeros((pathlength,T.shape[0]))
    Tori = T
    #M[0,:] = np.zeros((1,T.shape[0]))
    #M[0,start] = 1
    M[0,:] = start
    dispmatrix = M
    for i in range(1,pathlength):
        M[i,:] = np.dot(T.transpose(),M[0,:])
        dispmatrix[i,:] = M[i,:]*M[i,:]
        dispmatrix[i,:] = dispmatrix[i,:] / sum(dispmatrix[i,:])
        T = np.dot(T,Tori) 
    plt.imshow(dispmatrix.transpose())
    plt.show()
    return

    
def plot_TR_matrix(M, normalize = False):
# This method computes the transition rates  between all states of the Markov model.
# It assembles a matrix which in entry ij contains the TR of state i to state j. The matrix
# is then plotted. By setting the bool 'normalize' to 'True', the rows of the matrix are
# normalized
    n,m = M.T.shape
    TR = np.zeros((n,m))
    for i in range(0,n):
        for j in range(0,m):
            if (i is not j):
                disprobcur, effprobcur, FAB, KAB, TAB = M.tpt(np.array([i]),np.array([j]))
                TR[i,j] = KAB
    
    if (normalize):
        row_sums = TR.sum(axis=1)
        TR = TR / row_sums[:, np.newaxis]
    
    plt.imshow(TR)
    plt.show()
    return TR

def plot_stationary(stat):
    plt.bar(range(0,len(stat)),stat, alpha=0.5)
    plt.xlabel("State")
    plt.ylabel("Probability")
    plt.title("Stationary distribution")
    plt.show()

def plot_timescales(timescales):
    plt.plot(range(1,len(timescales)),timescales[1:],"*--")
    plt.xlabel("Eigenvalue $\mu_i$")
    plt.ylabel("Implied timescale")
    plt.show()
