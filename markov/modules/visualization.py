import matplotlib.pyplot as plt
import numpy as np

def plot_chain_frequencies(chain):
    states = list(set(chain))
    numofstates = len(states)
    count = np.zeros([numofstates,1])
    counter = 0
    
    for i in range(0,numofstates-1):
        count[counter] = (chain == states[i]).sum()
        counter = counter + 1
  
    plt.plot(sorted(count, reverse=True))
    plt.ylabel("Number of occurrences")
    plt.axis([0, max(states), 0, max(chain)])
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
