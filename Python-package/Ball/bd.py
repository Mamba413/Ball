# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 10:21:13 2019

@author: zhangyanhang
"""

from collections import namedtuple
from Ball.wrap_c import *
from Ball.utilize import *

bd_testResult = namedtuple('bd_testResult', ('statistic', 'pvalue'))
bd_Result = namedtuple('bd_Result', ('statistic'))

def bd_test(*args, **kwargs):
    """
    
    bd_test are ball divergence based multivariate nonparametric tests of two-sample
    or K-sample problem.
   
    Parameters
    ----------
    sample1, sample2, ... : array_like
        Two or more arrays with the sample measurements can be given asarguments.
    num_permutations : integer, optional
        the number of replications, when value equals to 0, the function returns the
        sample version of ball divergence. Default: num_permutations = 99.
    distance : 'True' or 'False', optional
        if distance = True, the input array will be considered as distance matrix.
        Default: distance = False.
    size : array_like
        a vector record sample size of each group. Needed only when the number of
        sample array is 1.
    weight : {'sum', 'max, 'maxsum'}, optional
        a character string used to choose the weight form of ball divergence statistic.
        Input must be one of the set above. Default: weight = sum
    
    Returns
    -------
    statistic : float
       ball divergence statistic
    pvalue : float
       The p-value for the Ball Divergence based permutation test

    Notes
    -----
    If only x is given, the statistic is computed from the original pooled samples, 
    stacked in matrix where each row is a multivariate observation, or from the distance matrix 
    when distance = TRUE. The first sizes[0] rows of x are the first sample, the next 
    sizes[1] rows of x are the second sample, etc.
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

    Test with ball divergence, when parameters x and y are one-dimensional array, and output its statistic and pvalue.

    >>> x = np.random.normal(0, 1, 50)
    >>> y = np.random.normal(1, 1, 50)
    >>> bd_test(x, y)
    bd_testResult(statistic=[0.196408479], pvalue=0.01)

    Test with ball divergence, when parameters x and y are multi-dimensional array, and output its statistic and pvalue.

    >>> x = np.random.normal(0, 1, 100).reshape(50, 2)
    >>> y = np.random.normal(3, 1, 100).reshape(50, 2)
    >>> bd_test(x, y)
    bd_testResult(statistic=[0.56810752], pvalue=0.01)

    Parameter size could input as array in ball divergence test, only and only if parameter y is not given.

    >>> n = 90
    >>> x = np.random.normal(0, 1, n)
    >>> bd_test(x, size=np.array([40, 50]))
    bd_testResult(statistic=[0.01993022499], pvalue=0.35)

    Parameter x could input as list in ball divergence test.

    >>> x = [np.random.normal(0, 1, num) for num in [40, 50]]
    >>> bd_test(x)
    bd_testResult(statistic=[0.013273874999], pvalue=0.47)

    When parameter x is a distance matrix, ball divergence test could also work on it, with parameter dst=False.

    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> sigma = [[1, 0], [0, 1]]
    >>> x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
    >>> y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
    >>> x = np.row_stack((x, y))
    >>> dx = euclidean_distances(x, x)
    >>> data_size = [50, 50]
    >>> bd_test(dx, size=data_size, dst=True)
    bd_testResult(statistic=[0.059688799999], pvalue=0.01)

    """
    num_permutations = 99
    num_thread = 1
    distance = False
    size = []
    weight = "sum"
    for key, value in kwargs.items():
        if key.lower() == "num_permutations" or key.lower() == "permutations":
            num_permutations = value
        elif key == "num_thread" or key.lower() == "nthread" or key.lower() == "thread":
            num_thread = value               
        elif key.lower() == "size":
            size = value
        elif key.lower() == "distance" or key.lower() == "dst":
            distance = value
        elif key.lower() == "weight":
            weight = value
        else:
            raise ValueError("input arguments is invalid")
            
    examine_permutations_arguments(num_permutations)
    examine_thread_arguments(num_thread)
    if len(args) == 1 and np.alen(args[0]) != 1 and distance == False and len(size) == 0:
        num_groups = len(args[0])
        args = args[0]
        map(examine_None, args)
    elif len(size) > 1:
        num_groups = len(size)        
    else:
        num_groups = len(args)
        map(examine_None, args)
        
    args = list(map(np.array, args))    
    if num_groups > 1:
        # return the dimension of samples
        p = examine_dimension(args)
        if len(size) == 0:
            size = list(map(np.alen, args))
        if p > 1:
            x = np.vstack(args)
            x = get_vectorized_distance_matrix(x)
            distance = True
        else:
            x = np.hstack(args)
            distance = False
    else:
        x = args[0]
        size = list(map(int, size))
        if len(x.shape) == 1:
            examine_size_arguments(x, size)
            distance = False
        else:
            examine_distance_matrix(x)
            examine_size_arguments(x, size)
            x = [x[i][j] for i in range(x.shape[0]) for j in range(x.shape[0]) if i<j]
            distance = True
        
    bd_stat, bd_pvalue = bd_test_wrap_c(x, size, num_permutations, distance, num_thread)
    
    if num_groups > 2:
        weight = weight.lower()
        if weight == "sum":
            bd_stat = bd_stat[0]
            bd_pvalue = bd_pvalue[0]
        elif weight == "max":
            bd_stat = bd_stat[1]
            bd_pvalue = bd_pvalue[1] 
        elif weight == "maxsum":
            bd_stat = bd_stat[2]
            bd_pvalue = bd_pvalue[2]
        else:
            raise ValueError("weight arguments is invalid!")        
    return bd_testResult(bd_stat, bd_pvalue)
                 
                 
def bd(*args, **kwargs):
    """

    calculate sample version of ball divergence

    Parameters
    ----------
    See bd_test.__doc__

    Returns
    -------
    statistic : float
       ball divergence statistic

    References
    ----------
    .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball
    divergence: nonparametric two sample test, \emph{The Annals of Statistics},
    to appear

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(7654567)   # fix seed to get the same result

    Test with ball divergence, when parameters x and y are one-dimensional array, and output its statistic.

    >>> x = np.random.normal(0, 1, 50)
    >>> y = np.random.normal(1, 1, 50)
    >>> bd(x, y)
    bd_Result(statistic=[0.196408479])

    Test with ball divergence, when parameters x and y are multi-dimensional array, and output its statistic.

    >>> x = np.random.normal(0, 1, 100).reshape(50, 2)
    >>> y = np.random.normal(3, 1, 100).reshape(50, 2)
    >>> bd(x, y)
    bd_Result(statistic=[0.56810752])

    Parameter size could input as array in ball divergence, only and only if parameter y is not given.

    >>> n = 90
    >>> x = np.random.normal(0, 1, n)
    >>> bd(x, size=np.array([40, 50]))
    bd_Result(statistic=[0.01993022499])

    Parameter x could input as list when compute ball divergence.

    >>> x = [np.random.normal(0, 1, num) for num in [40, 50]]
    >>> bd(x)
    bd_Result(statistic=[0.013273874999])

    When parameter x is a distance matrix, ball divergence could also work on it, with parameter dst=False.

    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> sigma = [[1, 0], [0, 1]]
    >>> x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
    >>> y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
    >>> x = np.row_stack((x, y))
    >>> dx = euclidean_distances(x, x)
    >>> data_size = [50, 50]
    >>> bd(x=dx, size=data_size, dst=True)
    bd_Result(statistic=[0.059688799999])

    """    
    for key in kwargs.keys():
        if key.lower() == "num_permutations" or key.lower() == "permutations":
            kwargs[key] = 0
    bd_value = bd_test(*args, **kwargs)
    return bd_Result(bd_value[0])
