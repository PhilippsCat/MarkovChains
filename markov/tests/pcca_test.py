import numpy as np
from markov.modules.analysis import MSM 
import markov.modules.estimate_transitionmatrix as estimate 
from nose.tools import assert_true

def test_pcca():
    n=1000
    k=5
    m=4
    chain   = np.random.randint(0,k,n)
    count   = estimate.count_transitions(chain)
    T       = estimate.db_estimator(count)
    M = MSM(T)
    assert_true(M.pcca(m).shape == (k,m))
