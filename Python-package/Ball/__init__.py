#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/7/26
# @Author  : Yanhang Zhang
# @Site    : 
# @File    : __init__.py.py

from ._cball import lib, ffi

# Optionally, provide Pythonic wrappers:
# def bcor_test(x, y, reps):
#     arr_x = ffi.new("double[]", x)
#     arr_y = ffi.new("double[]", y)
#     out = ffi.new("double *")
#     lib.bcor_test(arr_x, arr_y, out, len(x), len(y), reps)
#     return out[0]

# from Ball.bd import bd_test, bd
# from Ball.bcov import bcov_test, bcov
# from Ball.bcorsis import bcorsis, bcor
