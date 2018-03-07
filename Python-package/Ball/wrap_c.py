#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/2/5 20:26
# @Author  : Jin Zhu
# @Mail    : zhuj37@mail2.sysu.edu.cn
# @File    : wrap_c.py


from numpy import alen as npalen

from Ball.cball import bd_stat, bd_test
from Ball.cball import doubleArray, intArray


# BE CAUTIOUS!!! Init doubleArray for data with N element like this: doubleArray(N)
# One more again:
# doubleArray(N)   (Allow)
# doubleArray(N-1)  (Impermissible)


def bd_value_wrap_c(xy, size, weight, dst):
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
    num = npalen(xy)
    xy_copy = doubleArray(num)
    weight_copy = intArray(1)
    dst_copy = intArray(1)
    bd = doubleArray(1)
    # change the original data to doubleArray type
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
        pass
    weight_copy[0] = int(weight)
    dst_copy[0] = int(dst)
    bd[0] = 0.0
    #
    k = npalen(size)
    # calculate ball divergence statistic
    if k == 2:
        n1 = intArray(1)
        n2 = intArray(1)
        n1[0] = int(size[0])
        n2[0] = int(size[1])
        bd_stat(bd, xy_copy, n1, n2, weight_copy, dst_copy)
    else:
        pass
    bd_value = float(bd[0])
    return bd_value


def bd_test_wrap_c(xy, size, R, weight, dst):
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
    num = npalen(xy)
    xy_copy = doubleArray(num)
    weight_copy = intArray(1)
    dst_copy = intArray(1)
    R_copy = intArray(1)
    bd = doubleArray(1)
    permuted_bd = doubleArray(R)
    # change the original data to doubleArray type
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
        pass
    weight_copy[0] = int(weight)
    dst_copy[0] = int(dst)
    R_copy[0] = int(R)
    bd[0] = 0.0

    K = npalen(size)
    if K == 2:
        p = intArray(1)
        n1 = intArray(1)
        n2 = intArray(1)
        p[0] = int(1)
        n1[0] = int(size[0])
        n2[0] = int(size[1])
        bd_test(bd, permuted_bd, xy_copy, n1, n2, p, dst_copy, R_copy, weight_copy)
        # N = n1 + n2
    else:
        # size = int(size)
        # N = int(np.sum(size))
        # res <- .C("kbd_test", bd, permuted_bd, xy, size, N, K, dst, R, weight)
        pass

    # convert doubleArray to list:
    permuted_bd_list = []
    for i in range(R):
        permuted_bd_list.append(permuted_bd[i])
        pass

    bd_list = [bd[0]]
    return bd_list, permuted_bd_list

