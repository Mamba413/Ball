#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/7/24 20:26
# @Author  : Yanhang Zhang
# @Mail    : zhangyh93@mail2.sysu.edu.cn
# @File    : utilize.py


import numpy as np
import math
from sklearn.metrics.pairwise import euclidean_distances
from Ball.wrap_c import bcor_test_wrap_c


def examine_None(x):
    if x is not None:
        if np.any([np.isnan(i) for i in x]):
            raise ValueError("There is nan in the input!")
    else:
        raise ValueError("The input is None!")


def get_vectorized_distance_matrix(x):
    if len(x.shape) == 2:
        Dx = euclidean_distances(x, x)
        Dx = [Dx[i][j] for i in range(len(x)) for j in range(len(x)) if i < j]
        return Dx

    else:
        Dx = []
        for i in range(len(x)):
            for j in range(len(x)):
                if j > i:
                    Dx.append(abs(x[j] - x[i]))
        return np.array(Dx, dtype=float)


def examine_dimension(x):
    lens = []
    for sample in x:
        if len(sample.shape) == 1:
            lens.append(1)
        else:
            lens.append(sample.shape[1])
    if len(set(lens)) > 1:
        raise ValueError("Groups with different dimension!")
    return lens[0]


def examine_permutations_arguments(num):
    if num is None or num < 0:
        raise ValueError("num_permutations should be a value >=0!")


def examine_thread_arguments(num):
    if num is None or num < 0:
        raise ValueError("num_thread should be a value >=0!")


def examine_distance_matrix(x):
    dim = x.shape
    if dim[0] != dim[1]:
        raise ValueError("Distance matrix is not square")
    for i in range(dim[0]):
        for j in range(dim[0]):
            if x[i][j] < 0:
                raise ValueError("All elements should be positive!")
            elif i == j:
                if x[i][j] != 0:
                    raise ValueError("Diagonal of distance matrix is non-zero!")
            elif i < j:
                if round(x[i][j], 5) != round(x[j][i], 5):
                    print(i, j)
                    raise ValueError("Distance matrix is not symmetric")


def examine_size_arguments(x, size):
    if len(size) == 0:
        raise ValueError("size can't be None!")
    if np.any([x <= 0 for x in size]) or len(size) == 1:
        raise ValueError("size should be a list that length longer than 1, with all elements positive integer!")
    x_row = np.alen(x)
    n = np.sum(size)
    if x_row != n:
        raise ValueError("size arguments is invalid!")


def select_bcor_stat2(bcor_stat, weight):
    if type(weight) != str:
        raise ValueError("weight arguments is invalid")
    else:
        weight = weight.lower()
        if weight == "constant":
            bcor_stat = bcor_stat[0]
        elif weight == "probability" or weight == "prob":
            bcor_stat = bcor_stat[1]
        elif weight == "chisquare" or weight == "chisq":
            bcor_stat = bcor_stat[2]
        else:
            raise ValueError("weight arguments is invalid")
    return bcor_stat

def select_bcor_stat(y, x, x_num, f_num, n, distance, num_thread, weight):
    bcor_stat = bcor_test_wrap_c(y, x, x_num, f_num, n, distance, num_thread)
    if type(weight) != str:
        raise ValueError("weight arguments is invalid")
    else:
        weight = weight.lower()
        if weight == "constant":
            bcor_stat = bcor_stat[0]
        elif weight == "probability" or weight == "prob":
            bcor_stat = bcor_stat[1]
        elif weight == "chisquare" or weight == "chisq":
            bcor_stat = bcor_stat[2]
        else:
            raise ValueError("weight arguments is invalid")
    return bcor_stat


def get_x_index(ind, x_num):
    index = []
    for i in ind:
        if i == 0:
            temp = 0
        else:
            temp = np.sum(x_num[0:i])
        for i in range(x_num[i]):
            index.append(temp + i)
    return index


def examine_method_arguments(method):
    if type(method) != str or method.lower() not in ['standard', 'interaction', 'lm', 'gam']:
        raise ValueError('method arguments is invalid')


def examine_x_num_arguments(x, x_num):
    x_num = np.array(x_num)
    if len(x_num) <= 0:
        raise ValueError('x_num arguments is needed！')
    elif (x_num <= 0).any() == True:
        raise ValueError('All elements of x_num should be positive')
    elif (len(x) / np.sum(x_num) - int(len(x) / np.sum(x_num))) > 0:
        raise ValueError('x_num arguments is invalid')


def examine_params_arguments(params):
    params = np.array(params)
    if (params <= 0).any():
        raise ValueError('All elements of params should be positive')


def select_d_arguments(n, d):
    if type(d) != str and type(d) != int:
        raise ValueError('The type of d arguments is invalid')
    elif type(d) == str:
        d = d.lower()
        if d not in ['small', 'large']:
            raise ValueError('d arguments is invalid')
        elif d == 'small':
            d = math.floor(n / math.log(n))
        else:
            d = n - 1
    elif type(d) == int and d <= 0:
        raise ValueError('d arguments is invalid')

    return d


def examine_x_y_arguments(x, y):
    if len(x.shape) != len(y.shape):
        raise ValueError("Groups with different dimension!")
    elif len(x.shape) == 2 and x.shape[1] != y.shape[1]:
        raise ValueError("Groups with different dimension!")
    elif x.shape[0] != y.shape[0]:
        raise ValueError("Groups with different sample size!")
    return x.shape[0]
