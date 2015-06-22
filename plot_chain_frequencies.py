def plot_chain_frequencies(chain):
    import matplotlib.pyplot as plt
    import numpy as np
    
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
    return count