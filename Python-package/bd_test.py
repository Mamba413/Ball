# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 16:47:14 2018

@author: 99493
"""

from collections import namedtuple
from utilize import *
from wrap_c import *

bd_testResult = namedtuple('bd_testResult', ('statistic', 'pvalue'))


def bd_test(x, y=None, r=99, dst=False, size=[]):
    """
    bd_test are ball divergence based multivariate nonparametric tests of two-sample or 
    K-sample problem.
   
    Parameters
    ----------
    x a list, numeric vector, matrix, data.frame. 
    y a list, numeric vector, matrix or data.frame.
    r the number of replications, when R equals to 0, the function returns the sample version of ball divergence. Default: R = 99.
    dst if dst = True, x will be considered as a distance matrix. Default: dst = False.
    size a vector record sample size of each group.
    seed the random seed (TODO).
    
    Returns
    -------
    statistic : float
       ball divergence statistic
    pvalue : float
       The p-value for the Ball Divergence based permutation test

    Notes
    -----
    If only x is given, the statistic is 
    computed from the original pooled samples, stacked in 
    matrix where each row is a multivariate observation, or from the distance matrix 
    when dst = TRUE. The first sizes[1] rows of x are the first sample, the next 
    sizes[2] rows of x are the second sample, etc.
    If x is a list, its elements are taken as the samples to be compared, 
    and hence have to be numeric data vectors, matrix or data.frame.

    References
    ----------
    .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball 
    divergence: nonparametric two sample test, \emph{The Annals of Statistics},
    to appear
    
    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(7654567)   # fix seed to get the same result

    Example 1:
    >>> x = np.random.normal(0, 1, 50)
    >>> y = np.random.normal(1, 1, 50)
    >>> bd_test(x=x, y=y)
    bd_testResult(statistic=[2.0], pvalue=0.02)

    Example 2:
    >>> x = np.random.normal(0, 1, 100).reshape(50, 2)
    >>> y = np.random.normal(3, 1, 100).reshape(50, 2)
    >>> bd_test(x=x, y=y)
    bd_testResult(statistic=[0.012760320000000097], pvalue=0.02)

    Example 3:
    >>> n = 90
    >>> bd_test(np.random.normal(0, 1, n), size=np.array([40, 50]))
    bd_testResult(statistic=[0.013917274999999955], pvalue=0.01)

    Example 4:
    >>> x = [np.random.normal(0, 1, num) for num in [40, 50]]
    >>> bd_test(x)

    Example 5:
    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> sigma = [[1, 0], [0, 1]]
    >>> x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
    >>> y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
    >>> x = np.row_stack((x, y))
    >>> dx = euclidean_distances(x, x)
    >>> data_size = [50, 50]
    >>> bd_test(x=dx, size=data_size, dst=True)
    bd_testResult(statistic=[0.11696415999999976], pvalue=0.01)

    """
    weight = False
    method = 'permute'
    R = r

    if x is None or y is None:

        # examine input arguments x and y:
        if type(x) is list:
            for i in range(0, len(x)):
                size.append(x[i].shape[0])
                x[i] = np.matrix(x[i])
                x[i] = np.transpose(x[i])
            a = x[0]
            for i in range(1, len(x)):
                a = np.row_stack((a, x[i]))
            x = a

        else:
            x = get_matrixed_x(x, y)

        # examine input arguments size:
        examine_size_arguments(x, size)

        if dst > 0:
            xy = np.array(x.flatten())
        else:
            x1 = np.matrix(x)
            p = len(x1[0])
            if p > 1:
                xy = np.vstack((x, y))
                xy = get_vectorized_distance_matrix(xy)
                dst = True
            else:
                xy = x
                xy = np.asarray(xy)
                if xy.ndim != 1:
                    xy = xy[:, 0]
                dst = False

    else:
        x = np.matrix(x)
        y = np.matrix(y)
        # examine dimension:
        p = examine_dimension(x, y)

        if p > 1:
            xy = np.vstack((x, y))
            xy = get_vectorized_distance_matrix(xy)
            dst = True
            n1 = x.shape[0]
            n2 = y.shape[0]
            size = np.array([n1, n2])
        else:
            n1 = x.shape[0]
            n2 = y.shape[0]
            xy = np.row_stack((x, y))
            xy = np.array(xy.flatten())[0, :]
            dst = False
            size = np.array([n1, n2])

    # memory protect step:
    # memoryAvailable(n = sum(size), funs = 'BD.test')

    # examine R arguments:
    if method == "approx":
        R = 0
    else:
        examine_R_arguments(R)

    # main:
    if R == 0:
        result = bd_value_wrap_c(xy, size, weight, dst)
        # approximately method:
        if method == "approx":
            # if result[["info"]][["K"]] == 2:
            #     pvalue = calculatePvalue(np.prod(size) * result[["statistic"]] / np.sum(size),
            #                      BDTestNullDistribution)
            # else :
            #     return(result[["statistic"]])
            pass
        else:
            return result
    else:
        # permutation method:
        # examine seed arguments:
        # set.seed(examine_seed_arguments(seed))
        # hypothesis test:
        bd_value, bd_permuted_value = bd_test_wrap_c(xy, size, R, weight, dst)
        pvalue = calculatePvalue(bd_value, bd_permuted_value)
        return bd_testResult(bd_value, pvalue)

