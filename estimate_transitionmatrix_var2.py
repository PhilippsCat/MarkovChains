import numpy as np
def count_transitions(chain):
    n_markov_states = max(chain) + 1
    count_matrix = np.zeros((n_markov_states,n_markov_states))   
    for i in xrange(1, chain.shape[0]):
        count_matrix[chain[i-1]][chain[i]] = count_matrix[chain[i-1]][chain[i]] + 1
    return count_matrix

def estimate_transitionmatrix(chain):
    count_matrix = count_transitions(chain)
    for i in range(0,count_matrix.shape[0]):
        if sum(count_matrix[i,:]) > 0:
            count_matrix[i,:] = np.true_divide(count_matrix[i,:],sum(count_matrix[i,:]))
    return count_matrix

def db_estimator(chain, epsilon=1e-8, max_iter=50):
    """ 
    Estimate a reversible transition matrix from a given chain
    """
    count_matrix = count_transitions(chain)
   # x_old = estimate_transitionmatrix(chain) 
    x_old = np.eye(max(chain) + 1)
    n = len(x_old)
    x_new = np.zeros([n,n])
    c_rows = np.sum(count_matrix,1)
    iterations = 0
    x_older = x_new
    while np.allclose(x_old,x_older,epsilon) == False:
        x_old_rows = np.sum(x_old,1)        
        for i in range(n):
            if x_old_rows[i] == 0:
                continue            
            for j in range(n):
                if (x_old_rows[j] == 0) or (c_rows[i] == 0 and c_rows[j] == 0):
                    continue                
                x_new[i,j] = (count_matrix[i,j]+count_matrix[j,i])/(c_rows[i]/x_old_rows[i]+c_rows[j]/x_old_rows[j])
        iterations = iterations + 1
        x_older = x_old
        x_old = x_new
        if iterations == max_iter:
            break
      
    for i in range(n):
        row_sum = np.sum(x_new[i,:])
        if row_sum != 0:
            x_new[i,:] = x_new[i,:] / row_sum
    #x_new = x_new / np.sum(x_new,1)[:, np.newaxis]
    return x_new
      

#    for iterate in range(0, max_iter):
#        if np.allclose(x_old,x_new, epsilon):   #abort if given tolerance is reached earlier
#            break
#        x_old_rows=np.sum(x_old,1)
#        for i in range(0,n):
#            for j in range(0,n):
#                # fixed point iteration
#                x_new[i,j] = (count_matrix[i,j]+count_matrix[j,i])/(c_rows[i]/x_old_rows[i]+c_rows[j]/x_old_rows[j])
#        x_old=x_new
#
#    x_new = x_new / np.sum(x_new,1)[:, np.newaxis]  # normalize
#    return x_new
