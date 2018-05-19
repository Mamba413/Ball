#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/3/24 16:14
# @Author  : Mamba
# @Site    : 
# @File    : bcov.py


from collections import namedtuple

from Ball.wrap_c import *

from Ball.utilize import *

bcov_testResult = namedtuple('bcov_testResult', ('statistic', 'pvalue'))
bcov_Result = namedtuple('bcov_Result', ('statistic',))
bcor_Result = namedtuple('bcor_Result', ('statistic',))


def bcov_test(x, y=None, R=99, dst=False, weight=False, seed=4, num_threads=2):
    """
        bcov_test are non-parametric tests of multivariate independence in Banach space.
        The test decision is obtained via permutation, with replicates.

        Parameters
        ----------
        x a numeric vector, matirx, data.frame or dist object or list contains numeric vector, matrix or data.frame.
        y a numeric vector, matirx, data.frame or dist object.
        r the number of replications, when R equals to 0, the function returns the sample version of ball divergence. Default: R = 99.
        dst if dst = True, x will be considered as a distance matrix. Default: dst = False.
        size a vector record sample size of each group.
        weight when weight = True, weighted ball covariance or weighted ball correlation is used instead of ball covariance
        or ball correlation. Default: weight = FALSE
        type If type = 'bcor', ball correlation will be used instead of ball covariance.(default type = 'bcov')
        method if method = 'permute', a permutation procedure will be carried out, if method = 'approx', the p-values based on approximate Ball Covariance distribution are given.(Test arguments)


        Returns
        -------
        statistic : float
          ball covariance or ball correlation statistic
        pvalue : float
           The p-value for the Ball Divergence based permutation test

        Notes
        -----
        If two samples are pass to arguments x and y, the sample sizes (i.e. number of rows or length of the vector)
        of the two variables must agree.  Moreover, data pass to x or y must not contain missing or infinite values.
        If we set dst = TRUE, arguments x, y can be a dist object or a symmetric numeric matrix recording distance between samples;
        otherwise, these arguments are treated as data.

        The bcov_test statistic is bcov or bcor which are dependence measure in Banach space.
        The bcor test statistic is based on the normalized coefficient of ball covariance.

        References
        ----------
        .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball
        divergence: nonparametric two sample test, \emph{The Annals of Statistics},
        to appear

        Examples
        --------
        >>>error = np.array([random.uniform(-0.3,0.3) for _ in range(50)])
        >>>x = np.array([random.uniform(0, 4*math.pi) for _ in range(50)])
        >>>y = math.cos(x) + error
        >>>bcov_test(x, y)

        """
    method = 'permute'
    type_copy = "bcov"
    result = bcov_test_internal_wrap(x=x, y=y, R=R, dst=dst, weight=weight, seed=seed, method=method,
                                     type_copy=type_copy, num_threads=num_threads)
    return bcov_testResult(result[0], result[1])


def bcov_test_internal(x, y=None, R=99, dst=False, weight=False, seed=4, method='permute', type_copy='bcov',
                       num_threads=1):
    x = np.asmatrix(x)
    y = np.asmatrix(y)
    if x.shape[0] == 1:
        x = x.T
    if y.shape[0] == 1:
        y = y.T
    x_y_info = examine_x_y(x, y)
    num = x_y_info[0]
    p = x_y_info[1]

    if dst is False:
        if p != 1:
            x = get_vectorized_distance_matrix(x)
            y = get_vectorized_distance_matrix(y)
            dst = True
    else:
        x = np.asarray(x)
        y = np.asarray(y)

    # memory protect step:
    # memoryAvailable(n=sum(size), funs='BD.test')

    # examine test type:
    type_copy = examine_type_arguments(type_copy)

    if R == 0:
        result = bcov_value_wrap_c(x=x, y=y, n=num, weight=weight, dst=dst, type_copy=type_copy)
        if method == "approx":
            # pvalue = calculatePvalue(result[["statistic"]] * result[["info"]][["N"]], BITestNullDistribution)
            pass
        else:
            return result
    else:
        # set.seed(seed = examine_seed_arguments(seed))
        if len(x.shape) != 1 or len(y.shape) != 1:
            x = np.asarray(x).flatten()
            y = np.asarray(y).flatten()
        bcov_value, bcov_permuted_value = bcov_test_wrap_c(x=x, y=y, n=num, R=R, weight=weight, dst=dst,
                                                           type_copy=type_copy, num_threads=num_threads)
        pvalue = calculatePvalue(bcov_value, bcov_permuted_value)
        return bcov_value, pvalue


def kbcov_test_internal(x, R=99, dst=False, weight=False, seed=4, method='permute', type_copy='bcov'):
    x = np.asmatrix(x)
    size_list = [x[i].shape[1] for i in range(x.shape[0])]
    num = np.unique(size_list)
    if len(num) > 1:
        raise ValueError("sample sizes in each list must equal!")

    if dst:
        pass
    else:
        a = [0 for i in range(0, x.shape[0])]
        for i in range(0, x.shape[0]):
            a[i] = euclidean_distances(x[i].T, x[i].T)
        x = np.array(a)

    var_num = x.shape[0]

    # compute statistic:
    num = num[0]
    stat_value = doubleArray(1)
    stat_value[0] = kbcov_stat(x=x, num=num, var_num=var_num, weight=weight, type_copy=type_copy)

    if R == 0:
        pass
    else:
        # permutation procedure:
        permuted_stat = doubleArray(R)
        for r in range(0, R):
            x_copy = x
            for v in range(1, var_num):
                index = random.sample([i for i in range(0, num)], num)
                x_copy[v] = x[v][index, index]
            permuted_stat[r] = kbcov_stat(x=x_copy, num=num, var_num=var_num, weight=weight, type_copy=type_copy)

        permuted_stat_list = []
        for i in range(R):
            permuted_stat_list.append(permuted_stat[i])
            pass

        stat_value_list = [stat_value[0]]
        # calculate pvalue:
        pvalue = calculatePvalue(stat_value_list, permuted_stat_list)

    return stat_value_list, pvalue


def kbcov_stat(x, num, var_num, weight, type_copy):
    value_diff = doubleArray(1)
    stat_value = doubleArray(1)
    value_diff[0] = 0
    stat_value[0] = 0

    for i in range(0, num):
        for j in range(0, num):
            # compute ball statistic for ball with sample i as center and radius is d(x_{i}, x_{j})
            # d(x_{i}, x_{j}) are stored in x[[v]][i, j], x = 1, ..., var_num
            all_in_ball_vec = np.array([1 for i in range(0, num)])
            compare_list = []
            for v in range(0, var_num):
                a = np.array(x[v][:, i] <= x[v][i, j])
                compare_list.append(a)
                all_in_ball_vec = all_in_ball_vec * compare_list[v]

            prop_in_ball_vec = []
            for k in range(0, len(compare_list)):
                prop_in_ball_vec.append(np.mean(compare_list[k]))
            value_diff[0] = (np.mean(all_in_ball_vec) - np.prod(prop_in_ball_vec)) ** 2
            if weight:
                value_diff[0] = value_diff[0] / np.prod(prop_in_ball_vec)
            # aggrate statistic value:
            stat_value[0] = stat_value[0] + value_diff[0]

    return stat_value[0]


def bcov_test_internal_wrap(x, y, R, dst, seed, weight, method, type_copy, num_threads):
    if type(x) is list and y is None:
        result = kbcov_test_internal(x=x, R=R, dst=dst, weight=weight, seed=seed, method=method, type_copy=type_copy)
    else:
        result = bcov_test_internal(x=x, y=y, R=R, dst=dst, weight=weight, seed=seed, method=method,
                                    type_copy=type_copy, num_threads=num_threads)
    return result


def bcov(x, y, dst=False, weight=False, num_threads=2):
    """
    bcov computes ball covariance, which are multivariate measures of dependence in Banach space.

    Parameters
    ----------
    x a numeric vector, matirx, data.frame or dist object or list contains numeric vector, matrix or data.frame.
    y a numeric vector, matirx, data.frame or dist object.
    dst if dst = True, x will be considered as a distance matrix. Default: dst = False.
    weight when weight = True, weighted ball covariance or weighted ball correlation is used instead of ball covariance
    or ball correlation. Default: weight = FALSE

    Returns
    -------
    statistic : float
      ball covariance statistic

    Notes
    -----
    The sample sizes (number of rows or length of the vector) of the two variables must agree, and samples must not contain missing values.
    If we set \code{dst = TRUE}, arguments \code{x}, \code{y} can be a \code{dist} object or a
    symmetric numeric matrix recording distance between samples;
    otherwise, these arguments are treated as data.

    References
    ----------
    .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball
    divergence: nonparametric two sample test, \emph{The Annals of Statistics},
    to appear

    Examples
    --------
    >>>np.random.seed(5)
    >>>n = 50
    >>>x = np.random.normal(0, 1, n)
    >>>y = np.random.normal(0, 1, n)
    >>>bcov(x, y)

    """
    bcov_value = bcov_test_internal_wrap(x=x, y=y, R=0, dst=dst, seed=0, weight=weight, method="permute",
                                         type_copy="bcov",
                                         num_threads=num_threads)
    return bcov_Result(bcov_value)


def bcor(x, y, dst=False, weight=False, num_threads=2):
    """
    bcor computes ball correlation statistics, which are multivariate measures of dependence in Banach space.

    Parameters
    ----------
    x a numeric vector, matirx, data.frame or dist object or list contains numeric vector, matrix or data.frame.
    y a numeric vector, matirx, data.frame or dist object.
    dst if dst = True, x will be considered as a distance matrix. Default: dst = False.
    weight when weight = True, weighted ball covariance or weighted ball correlation is used instead of ball covariance
    or ball correlation. Default: weight = FALSE

    Returns
    -------
    statistic : float
      ball correlation statistic

    Notes
    -----
    The sample sizes (number of rows or length of the vector) of the two variables must agree, and samples must not contain missing values.
    If we set \code{dst = TRUE}, arguments \code{x}, \code{y} can be a \code{dist} object or a
    symmetric numeric matrix recording distance between samples;
    otherwise, these arguments are treated as data.

    References
    ----------
    .. [1]Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball
    divergence: nonparametric two sample test, \emph{The Annals of Statistics},
    to appear

    Examples
    --------
    >>>np.random.seed(5)
    >>>n = 50
    >>>x = np.random.normal(0, 1, n)
    >>>y = np.random.normal(0, 1, n)
    >>>bcor(x, y)

    """
    bcor_value = bcov_test_internal_wrap(x=x, y=y, R=0, dst=dst, seed=0, weight=weight, method="permute",
                                         type_copy="bcor",
                                         num_threads=num_threads)
    return bcor_Result(bcor_value)
