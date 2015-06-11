import numpy as np
import numpy.linalg as linalg
class MSM:
    
    def __init__(self,T):
        self.T = T

    """ Return the stationary distribution of the given Markov chain.
        Therefore, we compute the eigenvector of the transposed transition
        matrix corresponding to eigenvalue one and normalize.
    """    
    def stationarydistribution(self):
        EW,EV = linalg.eig(self.T.transpose())
        # By theory, there is only one biggest eigenvalue 1. As the linalg.eig
        # function may not order the eigenvalues, we use argmax() to identify
        # the index of the corresponding eigenvector and divide by its sum
        # to have a probability distribution. 
        statdist = EV[:,EW.argmax()] / EV[:, EW.argmax()].sum()
        return statdist

    """ Return the eigenvalues of the transition matrix as a numpy.ndarray
    """ 
    def eigenvalues(self):
        return linalg.eigvals(self.T) 
    
    # IS THiS REALLY WHAT WE WANT? Suppose only 2nd time scale is of interest
    def timescales(self):
        EW = self.eigenvalues()
        tscale = -1*np.ones(len(EW))/np.log(EW)
        return tscale

    def TPT(self, set):
        qminus = np.ones(len(T))
        qplus  = np.ones(len(T))
        n = len(qminus)
        f = np.zeros(shape=(n,n)) # discrete prob. curr.
        pi = self.stationarydistribution()
        effprobcur = np.zeros(shape=(n,n))
        FAB = 0
        for i in range(0,n):
            for j in range(0,n):
                if (i != j):

                    f[i,j] = pi[i]*qminus[i]*T[i,j]*qplus[j]
                    f[j,i] = pi[j]*qminus[j]*T[j,i]*qplus[i]
                    effprobcur[i,j] = max(0,(f[i,j]- f[j,i]))

                    if (i in set):
                        FAB = FAB + effprobcur[i,j]

        KAB = FAB / dot(pi,qminus)
        return f, effprobcur, FAB, KAB
    
    
    def hitting_prob(self, Set_A, Set_B = []):
        """ function gets a numpy.matrix and a hitting set Set_A [x,y,z...] as well as a 
        'Non hitting set' Set_B for computing comittors"""    
        coefficients = []
        numberOfStates = len(self.T)
        results = [0] * numberOfStates
        
        #take over all transition probabilities
        for i in range(numberOfStates):
            a = self.T[i]
            a[i] = a[i] - 1
            coefficients.append(a)  
            
        #set the probabilities for states in A to 1
        for i in range(len(Set_A)):
            b = [0] * numberOfStates
            b[Set_A[i]] = 1
            coefficients[Set_A[i]] = b
            results[Set_A[i]] = 1
        
        #set the probabilities for states in B to 0    
        for i in range(len(Set_B)):
            c = [0] * numberOfStates
            c[Set_B[i]] = 1
            coefficients[Set_B[i]] = c
     
        #make sure that the probabilities for absorbing nodes is 0 (as long as it is not in A) 
        for i in range(len(coefficients)):
            if coefficients[i] == [0] * numberOfStates:
                coefficients[i][i] = 1

        hitting_prob = np.linalg.solve(coefficients, results)
     
        return hitting_prob