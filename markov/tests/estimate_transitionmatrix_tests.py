""" Tests properties of estimate_transitionmatrix.
"""

import markov.modules.estimate_transitionmatrix as estimate
from nose.tools import assert_true 
import numpy as np

chain1 = np.array([0,1,1,0,2,0,1,0,2])
countMat1 = np.array([[0,2,2],
                      [2,1,0],
                      [1,0,0]])
transMat1 = np.array([[0,0.5,0.5],
                      [0.6666666666666666,0.3333333333333333,0],
                      [1,0,0]])

def count_transitions_test1():
    countMat = estimate.count_transitions(chain1)
    assert_true(np.array_equal(countMat,countMat1))

def estimate_transitionmatrix_test1():
    transMat = estimate.estimate_transitionmatrix(countMat1)
    assert_true(np.array_equal(transMat,transMat1))
