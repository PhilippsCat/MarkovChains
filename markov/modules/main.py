""" simple script that runs the pipline.
    Input: file path of data (as txt, whitespace seperated matrix)
"""
import markov.modules.clustering as cluster
import markov.modules.estimate_transitionmatrix as estimator
import markov.modules.analysis as analysis
import numpy as np
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--Ncc", dest="Ncc", default=2, type=int,
        help="number of cluster centers")
    p.add_argument("--dbMaxIter", dest="dbMaxIter", default=5000, type=int,
        help="maximum number of iterations for the fixpoint calculation of the db estimator")
    p.add_argument("--use-dbEst", dest="use_dbEst", default=False, action="store_true", 
        help="use db estimator to produce reversible matrix")
    p.add_argument("dataPath", 
        help="path to data matrix")
    args = p.parse_args()

    # run the pipeline (clustering, estimation, analysis)
    rawData = np.loadtxt(args.dataPath, ndmin=2)
    pt_to_cc, cc = cluster.kmeans(rawData, args.Ncc)
    countMat = estimator.count_transitions(pt_to_cc)
    if args.use_dbEst:
        transMat = estimator.db_estimator(countMat, max_iter=args.dbMaxIter)
    else:
        transMat = estimator.estimate_transitionmatrix(countMat)
    anaMat = analysis.MSM(transMat)

    print "stationary distribution: \n" + str(anaMat.statDist())
    print "eigenvalues: \n" + str(anaMat.eigenvalues())

    return True

if __name__ == '__main__':
    main()
