#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/2/5 20:26
# @Author  : Jin Zhu
# @Mail    : zhuj37@mail2.sysu.edu.cn
# @File    : wrap_c.py

from numpy import alen as npalen
from numpy import sum as npsum
from numpy import array as nparray
from numpy import column_stack
from numpy import double
from Ball.cball import bd_test, bcor_test, bcov_test, kbcov_test 
from Ball.cball import doubleArray, intArray



def bd_test_wrap_c(xy, size, R, distance, nthread=1):
    bd = doubleArray(6)
    pvalue = doubleArray(6)
    xy = nparray(xy, dtype=double)
    num = npalen(xy)
    
    xy_copy = doubleArray(num)
    num = npalen(size)
    size_copy = intArray(num)
    distance_copy = intArray(1)
    n = intArray(1)
    k = intArray(1)
    R_copy = intArray(1)
    nthread_copy = intArray(1)
    # change the original data to doubleArray type
    for i, xy_value in enumerate(xy):
        xy_copy[i] = xy_value
    for i, size_value in enumerate(size):
        size_copy[i] = int(size_value)
    n[0] = int(npsum(size))
    k[0] = int(npalen(size))
    distance_copy[0] = int(distance)
    R_copy[0] = int(R)
    nthread_copy[0] = int(nthread)
    # ball divergence based test
    bd_test(bd, pvalue, xy_copy, size_copy, n, k, distance_copy, R_copy, nthread_copy)
    # convert doubleArray to list:
    if k[0] > 2:
        pvalue_list = [pvalue[0], pvalue[2], pvalue[4]]
        bd_list = [bd[0], bd[2], bd[4]]
    else:
        pvalue_list = pvalue[0]
        bd_list = bd[0]
    return bd_list, pvalue_list


def bcor_test_wrap_c(y, x, x_num, f_num, n, dst_y, nthread):
    bcor_stat = doubleArray(3*f_num)
    y = nparray(y, dtype=double)
    x = nparray(x, dtype=double)
    
    y_copy = doubleArray(len(y))
    x_copy = doubleArray(len(x))
    x_num_copy = intArray(len(x_num))
    f_num_copy = intArray(1)
    size = intArray(1)
    n_copy = intArray(1)
    p = intArray(1)
    k = intArray(1)
    dst_y_copy = intArray(1)
    dst_x = intArray(1)
    nthread_copy = intArray(1)
    # change the original data to doubleArray type    
    for i, y_value in enumerate(y):
        y_copy[i] = y_value
        pass
    for i, x_value in enumerate(x):
        x_copy[i] = x_value
        pass
    for i,x_num_value in enumerate(x_num):
        x_num_copy[i] = x_num_value
    f_num_copy[0] = f_num
    size[0] = 0
    n_copy[0] = n
    p[0] = 0
    k[0] = 0
    dst_y_copy[0] = dst_y
    dst_x[0] = 0
    nthread_copy[0] = nthread
    
    bcor_test(bcor_stat, y_copy, x_copy, x_num_copy, f_num_copy, size, n_copy, p, k, dst_y_copy, dst_x, nthread_copy)
    # convert doubleArray to list:    
    bcor_stat_list = [bcor_stat[j] for j in range(3*f_num)]
    bcor_stat_list = column_stack(nparray(bcor_stat_list).reshape(f_num, 3))
    return bcor_stat_list


def bcov_test_wrap_c(x, y, n, R, distance, nthread):
    bcov_stat = doubleArray(3)
    pvalue = doubleArray(3)
    y = nparray(y, dtype=double)
    x = nparray(x, dtype=double)
    
    y_copy = doubleArray(len(y))
    x_copy = doubleArray(len(x))
    n_copy = intArray(1)
    R_copy = intArray(1)
    distance_copy = intArray(1)
    nthread_copy = intArray(1)
    # change the original data to doubleArray type
    for i, x_value in enumerate(x):
        x_copy[i] = x_value
    for i, y_value in enumerate(y):
        y_copy[i] = y_value
    n_copy[0] = int(n)
    distance_copy[0] = int(distance)
    R_copy[0] = int(R)
    nthread_copy[0] = int(nthread)
    
    bcov_test(bcov_stat, pvalue, x_copy, y_copy, n_copy, R_copy, distance_copy, nthread_copy)
    # convert doubleArray to list:    
    bcov_stat_list = [bcov_stat[0], bcov_stat[1], bcov_stat[2]]
    pvalue_list = [pvalue[0], pvalue[1], pvalue[2]]
    return bcov_stat_list, pvalue_list

def kbcov_test_wrap_c(x, k, n, R, nthread):
    kbcov_stat = doubleArray(3)
    pvalue = doubleArray(3)
    x = nparray(x, dtype=double)
    
    x_copy = doubleArray(len(x))
    k_copy = intArray(1)
    n_copy = intArray(1)
    R_copy = intArray(1)
    distance = intArray(1)
    nthread_copy = intArray(1)
    # change the original data to doubleArray type    
    for i, x_value in enumerate(x):
        x_copy[i] = x_value
    k_copy[0] = int(k)
    n_copy[0] = int(n)
    distance[0] = int(1)
    R_copy[0] = int(R)
    nthread_copy[0] = int(nthread)    
    
    kbcov_test(kbcov_stat, pvalue, x_copy, k_copy, n_copy, R_copy, distance, nthread_copy)
    # convert doubleArray to list:    
    kbcov_stat_list = [kbcov_stat[0], kbcov_stat[1], kbcov_stat[2]]
    pvalue_list = [pvalue[0], pvalue[1], pvalue[2]]
    return kbcov_stat_list, pvalue_list    