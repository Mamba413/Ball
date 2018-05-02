#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/2/5 20:26
# @Author  : Jin Zhu
# @Mail    : zhuj37@mail2.sysu.edu.cn
# @File    : wrap_c.py

from numpy import alen as npalen
from numpy import sum as npsum
from numpy import array as nparray
from numpy import double

from Ball.cball import bd_stat, bd_test
from Ball.cball import bcov_stat, bcov_test
from Ball.cball import doubleArray, intArray


# BE CAUTIOUS!!! Init doubleArray for data with N element like this: doubleArray(N)
# One more again:
# doubleArray(N)   (Allow)
# doubleArray(N-1)  (Impermissible)


def bd_value_wrap_c(xy, size, weight, dst, num_thread=1):
    '''
    wrapper C function which compute Ball Divergence

    Parameters
    ----------
    :param xy: a vectorized data (1D)
    :param size: a list object record sample size of each group.
    :param weight: A logical variable
    :param dst: A logical variable

    Returns
    -------
    bd : float
        Empirical Ball Divergence Value

    Examples
    ---
    >>> import numpy as np
    >>> x1 = np.random.normal(0, size=30)
    >>> x2 = np.random.normal(1, size=30)
    >>> x = np.append(x1, x2)
    >>> bd = bd_value_wrap_c(x, [30, 30], False, False)

    '''

    # initial C type structure
    xy = nparray(xy, dtype=double)
    num = npalen(xy)
    xy_copy = doubleArray(num)
    num = npalen(size)
    size_copy = intArray(num)
    weight_copy = intArray(1)
    dst_copy = intArray(1)
    bd = doubleArray(1)
    n = intArray(1)
    k = intArray(1)
    nthread = intArray(1)
    # change the original data to doubleArray type
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
        pass
    for i, size_value in enumerate(size):
        size_copy[i] = int(size_value)
        pass
    n[0] = int(npsum(size))
    k[0] = int(npalen(size))
    weight_copy[0] = int(weight)
    dst_copy[0] = int(dst)
    bd[0] = 0.0
    nthread[0] = int(num_thread)
    # calculate ball divergence statistics:
    bd_stat(bd, xy_copy, size_copy, n, k, weight_copy, dst_copy, nthread)
    bd_value = float(bd[0])
    return bd_value


def bd_test_wrap_c(xy, size, R, weight, dst, num_thread=1):
    '''

    Examples
    ---
    >>> import numpy as np
    >>> x1 = np.random.normal(0, size=30)
    >>> x2 = np.random.normal(1, size=30)
    >>> x = np.append(x1, x2)
    >>> bd_test_wrap_c(xy=x, size=[30, 30], R=10, weight=False, dst=False)

    '''

    # initial C type structure
    bd = doubleArray(1)
    permuted_bd = doubleArray(R)
    xy = nparray(xy, dtype=double)
    num = npalen(xy)
    xy_copy = doubleArray(num)
    num = npalen(size)
    size_copy = intArray(num)
    weight_copy = intArray(1)
    dst_copy = intArray(1)
    n = intArray(1)
    k = intArray(1)
    r = intArray(1)
    nthread = intArray(1)
    # change the original data to doubleArray type
    bd[0] = 0.0
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
        pass
    for i, size_value in enumerate(size):
        size_copy[i] = int(size_value)
        pass
    n[0] = int(npsum(size))
    k[0] = int(npalen(size))
    weight_copy[0] = int(weight)
    dst_copy[0] = int(dst)
    r[0] = int(R)
    nthread[0] = num_thread
    # ball divergence based test
    bd_test(bd, permuted_bd, xy_copy, size_copy, n, k, dst_copy, r, weight_copy, nthread)

    # convert doubleArray to list:
    permuted_bd_list = []
    for i in range(R):
        permuted_bd_list.append(permuted_bd[i])
        pass

    bd_list = [bd[0]]
    return bd_list, permuted_bd_list


def bcov_value_wrap_c(x, y, n, weight, dst, type_copy, num_threads=1):
    bcov = doubleArray(1)
    x = nparray(x, dtype=double)
    y = nparray(y, dtype=double)
    num = npalen(x)
    x_copy = doubleArray(num)
    for i, x_value in enumerate(x):
        x_copy[i] = x_value[0]
        pass
    num = npalen(y)
    y_copy = doubleArray(num)
    for i, y_value in enumerate(y):
        y_copy[i] = y_value[0]
        pass

    dst_copy = intArray(1)
    n_copy = intArray(1)
    weight_copy = intArray(1)
    num_threads_copy = intArray(1)

    dst_copy[0] = int(dst)
    n_copy[0] = int(n)
    weight_copy[0] = int(weight)
    num_threads_copy[0] = int(num_threads)

    if type_copy == "bcov":
        type1 = 1
    else:
        type1 = 2
    type_copy = intArray(1)
    type_copy[0] = int(type1)

    bcov_stat(bcov, x_copy, y_copy, n_copy, weight_copy, dst_copy, type_copy, num_threads_copy)
    bcov_value = float(bcov[0])
    return bcov_value


def bcov_test_wrap_c(x, y, n, R, weight, dst, type_copy, num_threads=1):
    bcov = doubleArray(1)
    permuted_bcov = doubleArray(R)
    x = nparray(x, dtype=double)
    y = nparray(y, dtype=double)
    num = npalen(x)
    x_copy = doubleArray(num)
    for i, x_value in enumerate(x):
        x_copy[i] = x_value
        pass
    num = npalen(y)
    y_copy = doubleArray(num)
    for i, y_value in enumerate(y):
        y_copy[i] = y_value
        pass

    n_copy = intArray(1)
    weight_copy = intArray(1)
    dst_copy = intArray(1)
    nthreads = intArray(1)
    r = intArray(1)

    n_copy[0] = int(n)
    weight_copy[0] = int(weight)
    dst_copy[0] = int(dst)
    nthreads[0] = int(num_threads)
    r[0] = int(R)
    if type_copy == "bcov":
        type1 = 1
    else:
        type1 = 2
    type_copy = intArray(1)
    type_copy[0] = int(type1)

    bcov_test(bcov, permuted_bcov, x_copy, y_copy, n_copy, r, weight_copy, dst_copy, type_copy, nthreads)

    permuted_bcov_list = []
    for i in range(R):
        permuted_bcov_list.append(permuted_bcov[i])
        pass

    bcov_list = [bcov[0]]
    return bcov_list, permuted_bcov_list
