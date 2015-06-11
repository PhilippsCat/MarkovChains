""" Tests for the markov property, irreducibility, periodicity and reversibility
    of a matrix or the corresponding MSM object.
""" 

import markov.modules.analysis as ana
from nose.tools import raises
import numpy as np

@raises(ana.InvalTransMatError)
def test_Reducibility1():
    ana.MSM(np.array([[1,0],[0,1]]))

@raises(ana.InvalTransMatError)
def test_Periodicity1():
    ana.MSM(np.array([[0,1],[1,0]]))

@raises(ana.InvalTransMatError)
def test_MarkovProperty1():
    ana.MSM(np.array([[1,0],[1,0]]))
