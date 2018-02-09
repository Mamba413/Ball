#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/2/5 20:26
# @Author  : Jin Zhu
# @Mail    : zhuj37@mail2.sysu.edu.cn
# @File    : wrap_c.py


from numpy import alen as npalen

from Ball.cball import bd_stat, bd_test
from Ball.cball import doubleArray


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
    import numpy as np
    x1 = np.random.normal(0, size=30)
    x2 = np.random.normal(1, size=30)
    x = np.append(x1, x2)
    bd = bd_value_wrap_c(x, [30, 30], False, False)

    '''

    # change the original data to doubleArray type
    num = npalen(xy)
    xy_copy = doubleArray(num)
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
        pass
    #
    bd = 0.0
    weight = int(weight)
    dst = int(dst)
    k = npalen(size)
    # calculate ball divergence statistic
    if k == 2:
        n1 = int(size[0])
        n2 = int(size[1])
        bd = bd_stat(xy_copy, n1, n2, weight, dst)
    else:
        pass
    bd_value = float(bd)
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

    num = npalen(xy)
    xy_copy = doubleArray(num)
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
        pass
    dst = int(dst)
    R = int(R)
    weight = int(weight)

    bd = doubleArray(1)
    permuted_bd = doubleArray(R)
    K = npalen(size)
    if K == 2:
        p = int(1)
        n1 = int(size[0])
        n2 = int(size[1])
        bd_test(bd, permuted_bd, xy_copy, n1, n2, p, dst, R, weight)
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

