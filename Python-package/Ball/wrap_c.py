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

