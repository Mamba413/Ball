# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 16:47:14 2018

@author: 99493
"""

from collections import namedtuple

from Ball.wrap_c import *

from Ball.utilize import *

bd_testResult = namedtuple('bd_testResult', ('statistic', 'pvalue'))
bd_Result = namedtuple('bd_Result', ('statistic', ))


def bd_test(x, y=None, R=99, dst=False, size=[], seed=4, num_threads=2):
    """
    bd_test are ball divergence based multivariate nonparametric tests of two-sample or 
    K-sample problem.
   
    Parameters
    ----------
    x a list, numeric vector, matrix.
    y a list, numeric vector, matrix.
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
    If only x is given, the statistic is computed from the original pooled samples, stacked in
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

    Test with ball divergence, when parameters x and y are one-dimensional array, and output its statistic and pvalue.

    >>> x = np.random.normal(0, 1, 50)
    >>> y = np.random.normal(1, 1, 50)
    >>> bd_test(x=x, y=y)
    bd_testResult(statistic=[0.196408479999], pvalue=0.01)

    Test with ball divergence, when parameters x and y are multi-dimensional array, and output its statistic and pvalue.

    >>> x = np.random.normal(0, 1, 100).reshape(50, 2)
    >>> y = np.random.normal(3, 1, 100).reshape(50, 2)
    >>> bd_test(x=x, y=y)
    bd_testResult(statistic=[0.75], pvalue=0.32)

    Parameter size could input as array in ball divergence test, only and only if parameter y is not given.

    >>> n = 90
    >>> x = np.random.normal(0, 1, n)
    >>> bd_test(x, size=np.array([40, 50]))
    bd_testResult(statistic=[0.0199302249999], pvalue=0.27)

    Parameter x could input as list in ball divergence test.

    >>> x = [np.random.normal(0, 1, num) for num in [40, 50]]
    >>> bd_test(x)
    bd_testResult(statistic=[0.013273874999], pvalue=0.51)

    When parameter x is a distance matrix, ball divergence test could also work on it, with parameter dst=False.

    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> sigma = [[1, 0], [0, 1]]
    >>> x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
    >>> y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
    >>> x = np.row_stack((x, y))
    >>> dx = euclidean_distances(x, x)
    >>> data_size = [50, 50]
    >>> bd_test(x=dx, size=data_size, dst=True)
    bd_testResult(statistic=[0.059310879999], pvalue=0.01)

    """
    # x = x.copy()
    # size = size.copy()
    weight = False
    method = 'permute'


    # if type(size) == np.ndarray:
    #     size = size.tolist()

    if x is None or y is None:
        # examine input arguments x and y:
        # if type(x) == pd.core.frame.DataFrame:
        #     x = np.array(x)
        #     x = np.transpose(x)
        #     if len(size) != len(x):
        #         for i in range(0, x.shape[0]):
        #             size.append(x[i].shape[0])
        if type(x) == list:
            size1=[]
            for i in range(0, len(x)):
                x[i] = np.matrix(x[i], dtype=np.double)
                x[i] = x[i].T
                if len(size) != len(x):
                    size1.append(x[i].shape[0])
            size = size1
            a = x[0]
            for i in range(1, len(x)):
                a = np.row_stack((a, x[i]))
            x = a
            examine_None(x, y)
        else:
            examine_None(x, y)
            x = get_matrixed_x(x, y)



        # examine input arguments size:
        examine_size_arguments(x, size)
        size = np.array(size)

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
                if len(np.array(xy.flatten()).shape) == 1:
                    xy = np.array(xy.flatten())
                else:
                    xy = np.array(xy.flatten())[0,:]


    else:
        # if type(x) == pd.core.frame.DataFrame:
        #     x = np.array(x)
        #     x = np.transpose(x)
        # if type(y) == pd.core.frame.DataFrame:
        #     y = np.array(y)
        #     y = np.transpose(y)
        if type(x) == list:
            x = np.array(x, dtype=np.double)
        if type(y) == list:
            y = np.array(y, dtype=np.double)
        examine_None(x, y)

        x = np.matrix(x).T
        y = np.matrix(y).T
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


def bd(x, y=None, dst=False, size=[]):
    '''

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
    >>> bd(x=x, y=y)
    bd_Result(statistic=[0.196408479999])

    Test with ball divergence, when parameters x and y are multi-dimensional array, and output its statistic.

    >>> x = np.random.normal(0, 1, 100).reshape(50, 2)
    >>> y = np.random.normal(3, 1, 100).reshape(50, 2)
    >>> bd(x=x, y=y)
    bd_Result(statistic=[0.75])

    Parameter size could input as array in ball divergence, only and only if parameter y is not given.

    >>> n = 90
    >>> x = np.random.normal(0, 1, n)
    >>> bd(x, size=np.array([40, 50]))
    bd_Result(statistic=[0.019930224999])

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
    bd_Result(statistic=[0.059310879999])

    '''

    bd_value = bd_test(x=x, y=y, R=0, dst=dst, size=size)
    return bd_Result(bd_value)
