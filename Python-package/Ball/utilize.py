#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/2/5 20:26
# @Author  : Yue Hu & Jin Zhu
# @Mail    : zhuj37@mail2.sysu.edu.cn &
# @File    : wrap_c.py


import numpy as np
import random
import math
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances


def get_matrixed_x(x, y):
    if x is None and y is None:
        raise ValueError("x and y are all null!")
    if x is None:
        x = y
    return np.array(x)

def examine_None(x, y):
    if x is not None:
        x = np.array(x.flatten())
        if np.any([np.isnan(i) for i in x]):
            raise ValueError("There is None in x!")
    if y is not None:
        y = np.array(y.flatten())
        if np.any([np.isnan(i) for i in y]):
            raise ValueError("There is None in y!")


def examine_size_arguments(x, size):
    # self examine:
    if size is None:
        raise ValueError("Size can't be None!")

    b = []
    for i in range(len(size)):
        if type(size[i]) is int:
            b.append(int(size[i]))
        else:
            b.append(0)
        size[i] = b[i]

    if np.any([x <= 0 for x in size]) or len(size) == 1:
        raise ValueError("Size should be a list that length longer than 1, with all elements positive integer!")

    # examine the consistency between x and size:
    if len(x.shape) == 1:
        x_row = x.shape[0]
    else:
        x_row = 1
        for i in range(len(x.shape)):
            x_row = x_row*x.shape[i]

    n = np.sum(size)
    if x_row != n:
        raise ValueError("Size arguments is invalid!")


def get_vectorized_distance_matrix(x):
    Dxy = euclidean_distances(x, x)
    #Dxy = np.array(Dxy.flatten())[0]
    Dxy = np.array(Dxy.flatten())
    return Dxy


def examine_dimension(x, y):
    p1 = len(x[0])
    p2 = len(y[0])
    if p1 != p2:
        raise ValueError("x and y with different dimension!")
    return p1


def examine_R_arguments(R):
    if R is None or R < 0:
        raise ValueError("r should be a value >=0 !")


def calculatePvalue(statValue, NullDistribution):
    return (np.sum(statValue < NullDistribution) + 1) / (len(NullDistribution) + 1)


# def memoryAvailable(n = np.sum(size), funs = 'BD.test'):
#     sysname = Sys.info()[1]
#     # memory check only available for window platform
#     if sysname == "Windows" :
#         MemoryAvailable = True
#         sizeLevel = (n / 1000)**2
#     if  funs == 'UBI.test'  :
#         queryMemory = 0.05 * sizeLevel
#     elif  funs == 'BI.test' :
#         queryMemory < - 0.16 * sizeLevel
#     elif  funs == 'UBD.test' :
#         queryMemory < - 0.02 * sizeLevel
#     elif  funs == 'BD.test' :
#         queryMemory < - 0.08 * sizeLevel
#
#     sysMemory = memory.limit() / 1024
#     if queryMemory > sysMemory :
#     # MemoryAvailable <- FALSE
#         raise ValueError("Sample size too large and system memory is not available!")
#     else :
#         if n > 8000 :
#             raise Warning("You may suffer from memory insufficient!")

def examine_x_y(x, y):
    dim_x = [x.shape[0], x.shape[1]]
    dim_y = [y.shape[0], y.shape[1]]

    if dim_x[0] == 1 and dim_x[1] == 1:
        if np.array(x).flatten()[0] is None:
            raise ValueError("x is None!")
    if dim_y[0] == 1 and dim_y[1] == 1:
        if np.array(y).flatten()[0] is None:
            raise ValueError("y is None!")

    n = dim_x[0]

    if n != dim_y[0]:
        raise ValueError("x and y have different sample sizes!")
    a = np.array(y).flatten()
    for i in range(0, len(a)):
        a[i] = (a[i]==None)
    if np.any(a):
        raise ValueError("Missing data in y!")
    a = np.array(x).flatten()
    for i in range(0, len(a)):
        a[i] = (a[i] == None)
    if np.any(a):
        raise ValueError("Missing data in x!")
    if dim_x[1] == 1 and dim_y[1] == 1:
        p = 1
    else :
        p = -1

    return np.array([n, p])

def examine_type_arguments(type_copy):
    if type_copy is not 'bcov' and type_copy is not 'bcor':
        type_copy = 'bcov'
    return type_copy

# def examine_seed_arguments(seed):
#   if seed is None:
#     seed = runif(1 , 0, .Machine$integer.max)???
#   return seed






