import numpy as np
def count_transitions(chain):
    n_markov_states = len(set(chain)) 
    count_matrix = np.zeros((n_markov_states,n_markov_states))
    for i in xrange(1, chain.shape[0]):
        count_matrix[chain[i-1]][chain[i]] = count_matrix[chain[i-1]][chain[i]] + 1
    return count_matrix

def estimate_transitionmatrix(chain):
    count_matrix = count_transitions(chain)
    for i in range(0,count_matrix.shape[0]):
        count_matrix[i,:] = np.true_divide(count_matrix[i,:],sum(count_matrix[i,:]))
    return count_matrix
