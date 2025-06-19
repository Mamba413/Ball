# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:24:01 2019

@author:Yanhang zhang
"""

from sklearn import linear_model
from pygam import LinearGAM
import numpy as np
from .utilize import *


def bcorsis(y, x, x_num, d="small", weight="constant", method="standard", dst_y=False, params=(5, 5), num_thread=0):
    """
    
    Generic non-parametric sure independence screening (SIS) procedure based on Ball Correlation. 
    Ball correlation is a generic multivariate measure of dependence in Banach space.  
    
    Parameters
    ----------

    x : numeric matirx
        included n rows and p columns. Each row is an observation vector and each column corresponding
        to a explanatory variable, generally p >> n.
    y : numeric vector, matrix
    x_num : array_like
        the array or list contains all the dimension of each variable. The length of x_num is the number of the variables
    d : integer, optional
        the hard cutoff rule suggests selecting d variables. Setting d = "large" or d = "small" means n - 1
        or floor(n/log(n)) variables are selected. If d is a integer, d variables are selected. Default: d = "small".
    weight : {'constant', 'probability', 'chisquare'}, optional
        a character string used to choose the weight form of Ball Covariance statistic.. Input must be one of "constant", "probability", or "chisquare".
        Default: weight = FALSE.
    method : {'standard', 'interaction', 'lm', 'gam'}, optional
        specific method for the BCor-SIS procedure. It must be one of "standard", "interaction", "lm" or "gam". Setting method = "standard" means
        performing standard SIS procedure while the options "lm" and "gam" mean carrying out iterative SIS procedure with ordinary linear regression
        and generalized additive models, respectively. The options "interaction" is designed for detecting variables with potential linear interaction 
        and associated with left censored responses, respectively. Default: method = "standard".
    dst_y : 'True' or 'False', optional
        if dst_y = TRUE, y will be considered as a distance matrix.
        Default: distance = FALSE.
    params : array_like
        parameters list only available when method = "lm" or "gam".        
    num_threads : integer
        Number of threads. If num_threads = 0, then all of available cores will be used. Default num_threads = 0.
    
    Returns
    -------
    A list of target variables  

    Notes
    -----
    bcorsis simultaneously computing Ball Correlation statistics with "constant", "probability", and "chisquare" weights. 
    Users can get other Ball Correlation statistics with different weight.
    
    "gam" method is slower than other methods when the number of variables is large.
    
    References
    ----------
    .. [1]Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic Sure Independence Screening Procedure, Journal of the
    American Statistical Association   
    .. [2]Jin, Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2018). Ball: An R package for detecting distribution difference and
    association in metric spaces.

    Examples
    --------
    >>> np.random.seed(10)    
    >>> x = np.random.normal(0, 1, n*p).reshape((n, p))
    >>> error = np.random.normal(0, 1, n)
    >>> y = 3 * x[:, 1] * x[:, 5] * x[:, 10] + error
    >>> result = bcorsis(y, x, x_num, method = "interaction", d=10)
    [1, 10, 567, 5, 1661]
    
    >>> np.random.seed(1000)
    >>> n = 150
    >>> p = 3000
    >>> mean = np.zeros(p)
    >>> cov = np.array([0.5]*np.square(p)).reshape((p, p))
    >>> cov[np.diag_indices(3000)] = 1
    >>> x = np.random.multivariate_normal(mean, cov, n)
    >>> error = np.random.normal(0, 1, n)
    >>> y = 4*np.square(x[:, 2])+6*np.square(x[:, 1])+8*x[:, 3]-10*x[:,4]+error
    >>> x_num = np.ones(3000)
    >>> target = [4, 1, 924, 2, 692, 3, 400, 2241, 2839, 2194, 170]
    >>> result = bcorsis(y, x, x_num, method="lm", params = [5, 3], d = 11)
    [4, 1, 924, 2, 692, 3, 400, 2241, 2839, 2194, 170]
    
    """
    examine_None(x)
    examine_None(y)
    x = np.array(x).T.flatten()
    examine_x_num_arguments(x, x_num)
    examine_method_arguments(method)

    f_num = len(x_num)
    if type(x_num) == np.ndarray:
        x_num = x_num.astype(int).tolist()
    n = int(len(x) / (np.sum(x_num)))
    d = select_d_arguments(n, d)
    if len(y.shape) > 1 and dst_y == True:
        examine_distance_matrix(y)
        y = [y[i][j] for i in range(np.alen(y)) for j in range(np.alen(y)) if i < j]
    else:
        y = np.array(y).T.flatten()

    method = method.lower()
    if method == "standard":
        bcor_stat = select_bcor_stat(y, x, x_num, f_num, n, dst_y, num_thread, weight)
        ind = []
        for i in range(d):
            ind.append(np.argmax(bcor_stat))
            bcor_stat[ind[i]] = -1

    elif method == "interaction":
        x_square = np.square(x)
        bcor_stat = select_bcor_stat(y, x, x_num, f_num, n, dst_y, num_thread, weight)
        bcor_square_stat = select_bcor_stat(y, x_square, x_num, f_num, n, dst_y, num_thread, weight)
        ind1 = []
        ind2 = []
        for i in range(d):
            ind1.append(np.argmax(bcor_stat))
            ind2.append(np.argmax(bcor_square_stat))
            bcor_stat[ind1[i]] = -1
            bcor_square_stat[ind2[i]] = -1

        ind = [index for index in ind1 if index in ind2]

    else:
        examine_params_arguments(params)
        d1 = params[0]
        d2 = params[1]
        bcor_stat = select_bcor_stat(y, x, x_num, f_num, n, dst_y, num_thread, weight)
        ind_have = []
        for i in range(d1):
            ind_have.append(np.argmax(bcor_stat))
            bcor_stat[ind_have[i]] = -1
        ind_last = ind_have

        x = x.reshape((np.sum(x_num), n)).T
        y = y.reshape((int(len(y) / n), n)).T
        ind_rest = [i for i in range(f_num) if i not in ind_have]

        if method == "lm":
            regr = linear_model.LinearRegression()
            while len(ind_have) < d:
                ind1 = get_x_index(ind_have, x_num)
                ind2 = get_x_index(ind_last, x_num)
                ind3 = get_x_index(ind_rest, x_num)

                regr.fit(x[:, ind1], x[:, ind3])
                x_new = x[:, ind3] - regr.predict(x[:, ind1])
                regr.fit(x[:, ind2], y)
                y = y - regr.predict(x[:, ind2])

                x_new = np.array(x_new).T.flatten()
                y = np.array(y).T.flatten()
                temp_x_num = [x_num[i] for i in ind_rest]
                bcor_stat = select_bcor_stat(y, x_new, temp_x_num, len(ind_rest), n, dst_y, num_thread, weight)
                ind_last = []
                for i in range(d2):
                    temp_ind = np.argmax(bcor_stat)
                    ind_last.append(ind_rest[temp_ind])
                    bcor_stat[temp_ind] = -1
                ind_have = ind_have + ind_last
                ind_rest = [i for i in ind_rest if i not in ind_have]

            ind = ind_have
        else:
            while len(ind_have) < d:
                ind1 = get_x_index(ind_have, x_num)
                ind2 = get_x_index(ind_last, x_num)
                ind3 = get_x_index(ind_rest, x_num)

                gam = LinearGAM().fit(x[:, ind1], x[:, ind3[0]])
                x_new = gam.deviance_residuals(x[:, ind1], x[:, ind3[0]])
                for i in range(len(ind3)):
                    if i == 0: continue
                    gam = LinearGAM().fit(x[:, ind1], x[:, ind3[i]])
                    x_new = np.hstack((x_new, gam.deviance_residuals(x[:, ind1], x[:, ind3[i]])))

                gam = LinearGAM().fit(x[:, ind2], y[:, 0])
                y = gam.deviance_residuals(x[:, ind2], y[:, 0])
                for i in range(y.shape[1]):
                    if i == 0: continue
                    gam = LinearGAM().fit(x[:, ind2], y[:, i])
                    y = np.hstack((y, gam.deviance_residuals(x[:, ind2], y[:, i])))

                temp_x_num = [x_num[i] for i in ind_rest]
                bcor_stat = select_bcor_stat(y, x_new, temp_x_num, len(ind_rest), n, dst_y, num_thread, weight)
                ind_last = []
                for i in range(d2):
                    temp_ind = np.argmax(bcor_stat)
                    ind_last.append(ind_rest[temp_ind])
                    bcor_stat[temp_ind] = -1
                ind_have = ind_have + ind_last
                ind_rest = [i for i in ind_rest if i not in ind_have]

            ind = ind_have

    return ind


def bcor(x, y, distance=False, weight="constant"):
    """
    
    calculate sample version of ball covariance

    Parameters
    ----------
    See bcorsis.__doc__

    Returns
    -------
    A list of ball correlation.

    References
    ----------
    .. [1]Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic Sure Independence Screening Procedure, Journal of the
    American Statistical Association   
    .. [2]Jin, Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2018). Ball: An R package for detecting distribution difference and
    association in metric spaces.

    Examples
    --------     
    >>> from Ball import bcor
    >>> import numpy as np
    >>> np.random.seed(1234)    
    >>> num = 5
    >>> x = np.random.normal(0, 1, num)
    >>> y = np.random.normal(0, 1, num)
    >>> bcor(x, y)
    
    >>> ## distance matrix input
    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> x = np.array(x, ndmin=2).T
    >>> y = np.array(y, ndmin=2).T
    >>> x = euclidean_distances(x, x)
    >>> y = euclidean_distances(y, y)
    >>> bcor(x, y, distance=True)
    """
    examine_None(x)
    examine_None(y)
    x = np.array(x, ndmin=2)
    y = np.array(y, ndmin=2)
    p = max(np.size(x, 0), np.size(y, 0))
    n = np.size(x, 1)
    
    dst_y = intArray(1)
    dst_x = intArray(1)   
    dst_y[0] = 1
    dst_x[0] = 1
    if distance:
        index = np.triu_indices(p, 1)
        x = x[index]
        y = y[index]
    else:
        if p != 1:
            y = get_vectorized_distance_matrix(y)
            x = get_vectorized_distance_matrix(x)
        else:
            dst_y[0] = 0
            dst_x[0] = 0  
            x = x[0]
            y = y[0]
            pass
        pass
    
    bcor_stat = doubleArray(3)
    y_copy = doubleArray(len(y))
    x_copy = doubleArray(len(x))
    for i, y_value in enumerate(y):
        y_copy[i] = y_value
        pass
    for i, x_value in enumerate(x):
        x_copy[i] = x_value
        pass
    x_number = intArray(1)
    x_number[0] = 1
    f_number = intArray(1)
    f_number[0] = 1
    num = intArray(1)
    num[0] = n
    p_copy = intArray(1)
    p_copy[0] = p
    k = intArray(1)
    k[0] = 1
    nth = intArray(1)
    nth[0] = 1

    bcor_test(bcor_stat, y_copy, x_copy, x_number, f_number, 
              num, p_copy, k, dst_y, dst_x, nth)
    bcor_stat_list = [bcor_stat[j] for j in range(3)]
    bcor_stat = select_bcor_stat2(bcor_stat_list, weight)
    return bcor_stat
