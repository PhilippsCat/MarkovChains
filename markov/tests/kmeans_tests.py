""" Tests properties for kmeans clustering.
"""

import markov.modules.clustering as clustering 
from nose.tools import assert_true 
import numpy as np

def test_trivial_kmeans():
    """ Tests for n distinct data points, whether each data point in a n-means 
    clustering is associated to a different cluster
    """
    n=10
    data = np.arange(n)
    data.shape=(-1,1)
    clusters = clustering.kmeans(data, n)
    assert_true(np.allclose(np.sort(clusters[0]), np.arange(n)))
