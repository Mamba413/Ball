#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
High-level Python wrappers for Ball C functions via pybind11.
"""
import numpy as np
from Ball._cball import py_bd_test, py_bcov_test, py_kbcov_test, py_bcor_test


def bd_test_wrap_c(xy, size, R, distance, nthread=1):
    xy = np.ascontiguousarray(xy, dtype=np.float64)
    size = np.ascontiguousarray(size, dtype=np.int32)
    bd_arr, pv_arr = py_bd_test(xy, size, int(R), int(distance), int(nthread))
    return bd_arr.tolist(), pv_arr.tolist()


def bcor_test_wrap_c(y, x, x_num, f_num, n, dst_y, nthread=1):
    y = np.ascontiguousarray(y, dtype=np.float64)
    x = np.ascontiguousarray(x, dtype=np.float64)
    x_num = np.ascontiguousarray(x_num, dtype=np.int32)
    out = py_bcor_test(y, x, x_num,
                       int(f_num), int(n), 0, 0,
                       int(dst_y), 0, int(nthread), 0)
    return out.reshape(f_num, 3).T


def bcov_test_wrap_c(x, y, n, R, distance, nthread=1):
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    stat_arr, pv_arr = py_bcov_test(x, y, int(n), int(R), int(distance), int(nthread))
    return stat_arr.tolist(), pv_arr.tolist()


def kbcov_test_wrap_c(x, k, n, R, nthread=1):
    x = np.ascontiguousarray(x, dtype=np.float64)
    stat_arr, pv_arr = py_kbcov_test(x, int(k), int(n), int(R), 1, int(nthread))
    return stat_arr.tolist(), pv_arr.tolist()
