#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
High-level Python wrappers for Ball C functions via CFFI.
"""
import numpy as np
from Ball.cball import lib, ffi


def bd_test_wrap_c(xy, size, R, distance, nthread=1):
    """
    Wrapper for C function:
      void bd_test(double *bd, double *pvalue, double *xy, int *size,
                   int *n, int *k, int *distance, int *R, int *nthread);
    """
    xy = np.asarray(xy, dtype=np.double)
    size = np.asarray(size, dtype=np.int32)
    n = np.sum(size)
    k = size.size

    # Allocate C arrays
    bd_c     = ffi.new('double[?]', k * 3)
    pvalue_c = ffi.new('double[?]', k * 3)
    xy_c     = ffi.new('double[?]', xy.size)
    size_c   = ffi.new('int[?]', size.size)
    n_c      = ffi.new('int *', int(n))
    k_c      = ffi.new('int *', int(k))
    dist_c   = ffi.new('int *', int(distance))
    R_c      = ffi.new('int *', int(R))
    thr_c    = ffi.new('int *', int(nthread))

    # Copy data
    for i, v in enumerate(xy):
        xy_c[i] = v
    for i, v in enumerate(size):
        size_c[i] = int(v)

    # Call C function
    lib.bd_test(bd_c, pvalue_c, xy_c, size_c, n_c, k_c, dist_c, R_c, thr_c)

    # Extract results
    bd_list = [bd_c[i] for i in range(k * 3)]
    pvalue_list = [pvalue_c[i] for i in range(k * 3)]
    return bd_list, pvalue_list


def bcor_test_wrap_c(y, x, x_num, f_num, n, dst_y, nthread=1):
    """
    Wrapper for C function:
      void bcor_test(double *out, double *y, double *x, int *x_num,
                     int *f_num, int *n, int *p, int *k,
                     int *dst_y, int *dst_x, int *nthread);
    Returns a (f_num,3) array.
    """
    y = np.asarray(y, dtype=np.double)
    x = np.asarray(x, dtype=np.double)
    x_num = np.asarray(x_num, dtype=np.int32)

    # Allocate
    out_c      = ffi.new('double[?]', f_num * 3)
    y_c        = ffi.new('double[?]', y.size)
    x_c        = ffi.new('double[?]', x.size)
    xnum_c     = ffi.new('int[?]', x_num.size)
    fnum_c     = ffi.new('int *', int(f_num))
    n_c        = ffi.new('int *', int(n))
    p_c        = ffi.new('int *', 0)
    k_c        = ffi.new('int *', 0)
    dst_y_c    = ffi.new('int *', int(dst_y))
    dst_x_c    = ffi.new('int *', 0)
    thr_c      = ffi.new('int *', int(nthread))

    # Copy
    for i, v in enumerate(y): y_c[i] = v
    for i, v in enumerate(x): x_c[i] = v
    for i, v in enumerate(x_num): xnum_c[i] = int(v)

    # Call
    lib.bcor_test(out_c, y_c, x_c, xnum_c, fnum_c, n_c, p_c, k_c, dst_y_c, dst_x_c, thr_c)

    # Reshape
    arr = np.frombuffer(ffi.buffer(out_c), dtype=np.double)
    return arr.reshape(f_num, 3)


def bcov_test_wrap_c(x, y, n, R, distance, nthread=1):
    """
    Wrapper for C function:
      void bcov_test(double *out, double *pvalue, double *x, double *y,
                     int *n, int *R, int *distance, int *nthread);
    Returns (stat_list, pvalue_list).
    """
    x = np.asarray(x, dtype=np.double)
    y = np.asarray(y, dtype=np.double)

    out_c      = ffi.new('double[3]')
    pvalue_c   = ffi.new('double[3]')
    x_c        = ffi.new('double[?]', x.size)
    y_c        = ffi.new('double[?]', y.size)
    n_c        = ffi.new('int *', int(n))
    R_c        = ffi.new('int *', int(R))
    dist_c     = ffi.new('int *', int(distance))
    thr_c      = ffi.new('int *', int(nthread))

    for i, v in enumerate(x): x_c[i] = v
    for i, v in enumerate(y): y_c[i] = v

    lib.bcov_test(out_c, pvalue_c, x_c, y_c, n_c, R_c, dist_c, thr_c)

    stat_list = [out_c[i] for i in range(3)]
    pvalue_list = [pvalue_c[i] for i in range(3)]
    return stat_list, pvalue_list


def kbcov_test_wrap_c(x, k, n, R, nthread=1):
    """
    Wrapper for C function:
      void kbcov_test(double *out, double *pvalue, double *x,
                      int *k, int *n, int *R, int *distance, int *nthread);
    Returns (stat_list, pvalue_list).
    """
    x = np.asarray(x, dtype=np.double)

    out_c      = ffi.new('double[3]')
    pvalue_c   = ffi.new('double[3]')
    x_c        = ffi.new('double[?]', x.size)
    k_c        = ffi.new('int *', int(k))
    n_c        = ffi.new('int *', int(n))
    R_c        = ffi.new('int *', int(R))
    dist_c     = ffi.new('int *', 1)
    thr_c      = ffi.new('int *', int(nthread))

    for i, v in enumerate(x): x_c[i] = v

    lib.kbcov_test(out_c, pvalue_c, x_c, k_c, n_c, R_c, dist_c, thr_c)

    stat_list = [out_c[i] for i in range(3)]
    pvalue_list = [pvalue_c[i] for i in range(3)]
    return stat_list, pvalue_list
