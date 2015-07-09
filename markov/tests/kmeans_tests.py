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

def test_calculate_centroids1():
    """ compute one centroid in 2 dimensions
    """
    data = np.array([[0,0],[2,0],[0,2],[2,2]])
    means = np.array([[-1,-1]])
    point_to_means = np.zeros(4)
    center = clustering.calculate_centroids(data,means,point_to_means)
    assert_true(np.array_equal(center[0],[1,1]))

def test_closest_means1():
    """ compute the mapping with two clusters
        in 2 dimensions
    """
    data = np.array([[0,0],[2,0],[0,2],[2,2],
                     [7,0],[9,0],[7,2],[9,2]])
    means = np.array([[1,1],[8,1]])
    point_to_means = clustering.closest_means(data,means)
    assert_true(np.array_equal(point_to_means,[0,0,0,0,1,1,1,1]))
