#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/2/5 20:26
# @Author  : Yue Hu & Jin Zhu
# @Mail    : zhuj37@mail2.sysu.edu.cn &
# @File    : wrap_c.py


import numpy as np
from sklearn.metrics.pairwise import euclidean_distances


def get_matrixed_x(x, y):
    if x is None and y is None:
        raise ValueError("x and y are all null!")
    if x is None:
        x = y
    return np.array(x)


def examine_size_arguments(x, size):
    # self examine:
    if size is None:
        raise ValueError("size arguments is needed")

    size = [int(x) for x in size]
    if np.any(np.isnan(size)) or np.any([x <= 0 for x in size]) or len(size) == 1:
        raise ValueError("size arguments is invaild!")

    # examine the consistency between x and size:
    x_row = x.shape[0]
    n = np.sum(size)
    if x_row != n:
        raise ValueError("size arguments is invaild!")


def get_vectorized_distance_matrix(x):
    Dxy = euclidean_distances(x, x)
    Dxy = np.array(Dxy.flatten())[0, :]
    return Dxy


def examine_dimension(x, y):
    p1 = len(x[0])
    p2 = len(y[0])
    if p1 != p2:
        raise ValueError("x and y with different dimension!")
    return p1


def examine_R_arguments(R):
    if R is None or R < 0:
        raise ValueError("R arguments is invaild!")


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

