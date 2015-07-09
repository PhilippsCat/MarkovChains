import numpy as np
from nose.tools import assert_true
from markov.modules.clustering import *
def regspace_test():
    chain = np.array([0.9,1,1.1,1.9,2,2.1])
    regspace_states, regspace_centers = regspace(chain,0.5)
    print regspace_states
    assert_true(np.array_equal(regspace_states,[0.,0.,0.,1.,1.,1.]))
    
    
    
