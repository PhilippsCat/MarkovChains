import matplotlib.pyplot as plt
import numpy as np

#This method sorts the states by their number of occurrences in 
# a chain und plots the numbers of occurrences in descending order. It requires
# a 1D-numpy-array
def plot_chain_frequencies(chain):

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
