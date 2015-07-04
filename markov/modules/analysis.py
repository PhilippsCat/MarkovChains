import numpy as np
import numpy.linalg as linalg
import fractions as frac
from collections import deque

class InvalTransMatError(ValueError):
    """Raised in ctor of MSM when the matrix 
       doesn't conform to basic sanity checks.
    """
    pass

def markovProp(T):
    """Check if the matrix has the markov property.
    """
    if not (T.ndim == 2):
        return False, "wrong number of dimensions"
    d1, d2 = T.shape
    if not (d1 == d2):
        return False, "not a square matrix"
    if not np.allclose(np.sum(T,1), np.ones(len(T))):
        return False, "rows don't sum up to one"
    return True, "all good"

def commClasses(g):
    """ compute the communication classes
        of the graph, given as a transition
        matrix.
    """

    def dfs(graph, node, vis_nodes):
        """depth first search (non-destructive)
        """
        g = np.copy(graph)
        return dfs_destr(g, node, vis_nodes)


    def dfs_destr(graph, node, vis_nodes):
        """depth first search (destructive)
        """
        trans = graph[node] > 0
        # deleting inbound edges == marking as visited 
        graph[:, node] = 0
        ret = []
        for i in range(len(trans)):
            if (trans[i] 
                and (i not in vis_nodes)
                and (i != node)):
                    ret += dfs_destr(graph, i, vis_nodes)
        ret.append(node)
        return ret

    def check(V, n): 
        """checking for next free node in V
        """
        for node in range(n):
            if node not in V:
                return node
        return -1
    
    # avoiding destruction of graph
    graph = np.copy(g)
    communic_classes = []
    V = []
    n = len(graph)
    
    # explore all nodes
    node = check(V, n)
    while (node != -1):
        V += dfs(graph, node, V)
        node = check(V, n)
        
    # continue in transposed graph
    graph_t = graph.transpose()
    while V:
        node = V[-1]
        C = dfs(graph_t, node, [])
        communic_classes.append(C)
        for i in C:
            if i in V: V.remove(i)
            graph_t[i,:] = 0 
            graph_t[:,i] = 0
    
    return communic_classes

def period(g):
    """ computes the period of an irreducible 
        markov matrix.
        algorithm taken from:
        Graph-Theoretic Analysis of Finite Markov Chains,
        Jarvis and Shier,
        chapter 17
    """

    def val(a, b):
        return a - b + 1

    lvls = {0 : 0} # node 0 is on lvl 0
    Q = deque([0]) # queue of nodes to be processed
    per = 0        # no actual period is set yet 

    while Q:
        w = Q.pop()
        trans = g[w] > 0
        for v in range(len(trans)):             
            if trans[v]:
                if (v in lvls.keys()): # v has already been visited
                    p = val(lvls[w],lvls[v])
                    if (per == 0):
                        # set initial period
                        per = p
                    if (p == 1): 
                        return 1
                    if (p > 1): 
                        per = frac.gcd(per,p)
                else:
                    Q.appendleft(v)
                    lvls[v] = lvls[w] + 1

    return per

class MSM:
    
    def __init__(self,T):
        """In this implementation only irreducible, aperiodic
           and reversible markov matrices can be created.
           Trying to create a MSM object that does not fulfill
           those properties will lead to an InvalTransMatError.
        """
        mP, msg = markovProp(T)
        if not mP:
            raise InvalTransMatError("not a markov matrix: " + msg)
        self.T = T
        self.pi = self.statDist()
        if not (len(commClasses(T)) == 1):
            raise InvalTransMatError("matrix is reducible")
        if not (period(T) == 1):
            raise InvalTransMatError("matrix is aperiodic")
        PI = np.diag(self.pi)
        if not np.allclose(np.dot(PI, T), np.dot(T.transpose(),PI)): 
            # detailed balance condition
            raise InvalTransMatError("matrix is not reversible")
 
    def statDist(self):
        """ Return the stationary distribution of the given Markov chain.
            Therefore, we compute the eigenvector of the transposed transition
            matrix corresponding to eigenvalue one and normalize.
        """
        EW,EV = linalg.eig(self.T.transpose())
        # By theory, there is only one biggest eigenvalue 1. As the linalg.eig
        # function may not order the eigenvalues, we use argmax() to identify
        # the index of the corresponding eigenvector and divide by its sum
        # to have a probability distribution. 
        statdist = EV[:,EW.argmax()] / EV[:, EW.argmax()].sum()
        return statdist

    def eigenvalues(self):
        """ Return the eigenvalues of the transition matrix as a numpy.ndarray
        """
        return linalg.eigvals(self.T) 
    
    # IS THiS REALLY WHAT WE WANT? Suppose only 2nd time scale is of interest
    def timescales(self):
        EW = self.eigenvalues()
        tscale = -1*np.ones(len(EW))/np.log(EW)
        return tscale

    def tpt(self, setA, setB):
        n = self.T.shape[0]
        qminus = self.hitting_prob(setA, setB) # backward committor
        qplus  = self.hitting_prob(setB, setA) # forward committor
        f = np.zeros((n,n)) # discrete prob. curr.
        pi = self.pi # stationary distribution
        effprobcur = np.zeros((n,n))
        FAB = 0 # average total number of reactive trajectories
            for i in range(0,n):
                for j in range(0,n):
                    if (i != j):
                        f[i,j] = pi[i]*qminus[i]*T[i,j]*qplus[j]
                        f[j,i] = pi[j]*qminus[j]*T[j,i]*qplus[i]
                        effprobcur[i,j] = max(0,(f[i,j]- f[j,i]))

                        if (i in setA):
                            FAB = FAB + effprobcur[i,j]
            disprobcur = f  
        KAB = FAB / np.dot(pi,qminus) # transition rate
        TAB = 1/KAB #mean first passage time
        
        return disprobcur, effprobcur, FAB, KAB, TA
    
    
    def hitting_prob(self, Set_A, Set_B = []):
        """ function gets a numpy.matrix and a hitting set Set_A [x,y,z...] as well as a 
        'Non hitting set' Set_B for computing comittors
        it is assumed that the states starting with '1' and not '0'"""    
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
            b[Set_A[i] - 1] = 1
            coefficients[Set_A[i] - 1] = b
            results[Set_A[i] - 1] = 1
        
        #set the probabilities for states in B to 0    
        for i in range(len(Set_B)):
            c = [0] * numberOfStates
            c[Set_B[i] - 1] = 1
            coefficients[Set_B[i] - 1] = c
     
        #make sure that the probabilities for absorbing nodes is 0 (as long as it is not in A) 
        for i in range(len(coefficients)):
            if coefficients[i] == [0] * numberOfStates:
                coefficients[i][i] = 1

        hitting_prob = np.linalg.solve(coefficients, results)
     
        return hitting_prob
