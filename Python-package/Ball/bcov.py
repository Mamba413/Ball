# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:11:32 2019

@author: zhangyanhang
"""

from collections import namedtuple
from Ball.wrap_c import *
from Ball.utilize import *

bcov_testResult = namedtuple('bcov_testResult', ('statistic', 'pvalue'))
bcov_Result = namedtuple('bcov_Result', ('statistic'))

def bcov_test(*args, **kwargs):
    """
    bcov_test are non-parametric tests of multivariate independence in Banach space.
    The test decision is obtained via permutation, with replicates.

    Parameters
    ----------
    sample1, sample2, ... : array_like
        Two or more arrays with the sample measurements can be given asarguments.
    num_permutations : integer, optional
        the number of replications, when value equals to 0, the function returns the
        sample version of ball divergence. Default: num_permutations = 99.
    distance : 'True' or 'False', optional
        if distance = True, the input array will be considered as distance matrix.
        Default: distance = False.        size a vector record sample size of each group.
    weight : {'constant', 'probability', 'chisquare'}, optional
        a character string used to choose the weight form of ball covariance statistic.
        Input must be one of Input must be one of "constant", "probability", or "chisquare".
        Default: weight = constant

    Returns
    -------
    statistic : float
      ball covariance or ball correlation statistic
    pvalue : float
       The p-value for the Ball Divergence based permutation test

    Notes
    -----
    If k samples are pass to arguments args, the sample sizes (i.e. number of rows or length of the vector)
    of the k variables must agree.  Moreover, data must not contain missing or infinite values.
    If we set distance = TRUE, arguments args can be a dist object or a symmetric numeric matrix recording distance between samples;
    otherwise, these arguments are treated as data.

    The bcov_test statistic is bcov or bcor which are dependence measure in Banach space.

    References
    ----------
    .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball
    divergence: nonparametric two sample test, \emph{The Annals of Statistics},
    to appear

    Examples
    --------
    >>> np.random.seed(2)
    >>> error = np.random.uniform(-0.3, 0.3, 50)
    >>> x = np.random.uniform(0, 4*math.pi, 50)
    >>> y = np.cos(x) + error
    >>> bcov_test(x, y)
    bcov_testResult(statistic=[0.0021155349], pvalue=0.01)
    
    >>> x = np.random.normal(0, 1, 50)
    >>> y = [1 if i > 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
    >>> z = [1 if i < 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
    >>> bcov_test(x, y, z)
    bcov_testResult(statistic=[0.00252983469], pvalue=0.01)    
    
    """
    weight = "constant"
    num_permutations = 99
    num_thread = 2
    distance = False
    for key, value in kwargs.items():
        if key.lower() == "weight":
            weight = value
        elif key.lower() == "num_permutations" or key.lower() == "permutations":
            num_permutations = value
        elif key.lower() == "num_thread" or key.lower() == "nthread" or key.lower() == "thread":
            num_thread = value  
        else:
            raise ValueError("input arguments is invalid")            
    
    args = list(map(np.array, args))
    list(map(examine_None, args))
    num_groups = len(args) 
    n = list(map(np.alen, args))
    if len(set(n)) != 1:
        raise ValueError("Variables with different sample size")
    else:
        n = n[0]
        
    if num_groups == 2:
        dim = list(map(lambda x:len(x.shape), args))
        if dim[0] != 1 or dim[1] != 1:
            args = list(map(get_vectorized_distance_matrix, args))
            distance = True
        bcov_stat, bcov_pvalue = bcov_test_wrap_c(args[0], args[1], n, num_permutations, distance, num_thread)
    else:
        args = list(map(get_vectorized_distance_matrix, args))
        x = np.hstack(args)
        bcov_stat, bcov_pvalue = kbcov_test_wrap_c(x, num_groups, n, num_permutations, num_thread)
    
    if type(weight) != str:
        raise ValueError("weight arguments is invalid")
    else:   
        weight = weight.lower()
        if weight == "constant":
           bcov_stat = bcov_stat[0]
           bcov_pvalue = bcov_pvalue[0]
        elif weight == "probability" or weight == "prob":
           bcov_stat = bcov_stat[1]
           bcov_pvalue = bcov_pvalue[1]       
        elif weight == "chisquare" or weight == "chisq":
           bcov_stat = bcov_stat[2]
           bcov_pvalue = bcov_pvalue[2]
        else:
            raise ValueError("weight arguments is invalid!")

    return bcov_testResult(bcov_stat, bcov_pvalue)

def bcov(*args, **kwargs):
    """
    
    calculate sample version of ball covariance

    Parameters
    ----------
    See bcov_test.__doc__

    Returns
    -------
    statistic : float
      ball covariance or ball correlation statistic

    References
    ----------
    .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball
    divergence: nonparametric two sample test, \emph{The Annals of Statistics},
    to appear

    Examples
    --------
    >>> np.random.seed(2)
    >>> error = np.random.uniform(-0.3, 0.3, 50)
    >>> x = np.random.uniform(0, 4*math.pi, 50)
    >>> y = np.cos(x) + error
    >>> bcov(x, y)
    bcov_Result(statistic=[0.0021155349])
    
    >>> x = np.random.normal(0, 1, 50)
    >>> y = [1 if i > 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
    >>> z = [1 if i < 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
    >>> bcov(x, y, z)
    bcov_Result(statistic=[0.00252983469])      
    
    """
    for key in kwargs.keys():
        if key.lower() == "num_permutations" or key.lower() == "permutations":
            kwargs[key] = 0
    bcov_value = bcov_test(*args, **kwargs)
    return bcov_Result(bcov_value[0])