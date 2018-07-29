#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/3/8 16:36
# @Author  : Mamba
# @File    : test_bd.py


import unittest
from unittest import TestCase
from Ball.bd import bd_test, bd
import numpy as np
import pandas as pd


class Test_bd(TestCase):
    def test_bd(self):
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = np.array(x, dtype=np.double)
        bd_value = bd(x=x, size=np.array([5, 5]))
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.72319999999)

        np.random.seed(7654567)
        x = np.random.normal(0, 1, 50)
        y = np.random.normal(1, 1, 50)
        bd_value = bd(x=x, y=y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.196408479999)

        x = np.random.normal(0, 1, 100).reshape(50, 2)
        y = np.random.normal(3, 1, 100).reshape(50, 2)
        bd_value = bd(x=x, y=y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.75)

        n = 90
        x = np.random.normal(0, 1, n)
        bd_value = bd(x, size=np.array([40, 50]))
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.019930224999)

        x = [np.random.normal(0, 1, num) for num in [40, 50]]
        bd_value = bd(x=x)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.013273874999)

        from sklearn.metrics.pairwise import euclidean_distances
        sigma = [[1, 0], [0, 1]]
        x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
        y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
        x = np.row_stack((x, y))
        dx = euclidean_distances(x, x)
        data_size = [50, 50]
        bd_value = bd(x=dx, size=data_size, dst=True)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.059310879999)



    def test_bd_test(self):
        x = [1, 2, 3, 4, 5]
        y = [1, 2, 3, 4, 5]
        bd_value = bd_test(x, y)
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.0)

        np.random.seed(7654567)
        x = np.random.normal(0, 1, 50)
        y = np.random.normal(1, 1, 50)
        bd_value = bd_test(x=x, y=y)
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.196408479999)

        x = np.random.normal(0, 1, 100).reshape(50, 2)
        y = np.random.normal(3, 1, 100).reshape(50, 2)
        bd_value = bd_test(x=x, y=y)
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.75)

        n = 90
        x = np.random.normal(0, 1, n)
        bd_value = bd_test(x, size=np.array([40, 50]))
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.019930224999)

        x = [np.random.normal(0, 1, num) for num in [40, 50]]
        bd_value = bd_test(x=x)
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.013273874999)

        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = np.array(x, dtype=np.double)
        bd_value = bd_test(x=x, size=np.array([5, 5]))
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.7232)

        x = [1, 2, 3, 4, None, 5, 6, np.nan, 7, 9]
        self.assertRaises(ValueError, lambda: bd_test(x=x, size=np.array([5, 5])))

        x = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]
        bd_value = bd_test(x=x)
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 2.4032)

        x = [[1, 2, None, 4, 5], [6, 7, 8, np.nan, 10], [11, 12, 13, 14, 15]]
        self.assertRaises(ValueError, lambda: bd_test(x=x))

        # x = {"first": [1, 2, 3, 4, 5], "second": [6, 7, 8, 9, 10], "third": [11, 12, 13, 14, 15]}
        # x = pd.DataFrame(x)
        # bd_value = bd_test(x)
        # bd_value = bd_value[0][0]
        # self.assertAlmostEqual(bd_value, 2.4032)

        # x = {"first": [1, 2, None, 4, 5], "second": [6, 7, 8, np.nan, 10], "third": [None, 12, 13, 14, 15]}
        # x = pd.DataFrame(x)
        # self.assertRaises(ValueError, lambda: bd_test(x))
        #
        # x = {"first": [1, 2, 3, 4, 5], "second": [6, 7, 8, 9, 10]}
        # y = {"first": [11, 12, 13, 14, 15], "second": [16, 17, 18, 19, 20]}
        # x = pd.DataFrame(x)
        # y = pd.DataFrame(y)
        # bd_value = bd_test(x=x, y=y)
        # bd_value = bd_value[0][0]
        # self.assertAlmostEqual(bd_value, 0.875)
        #
        # x = {"first": [1, 2, None, 4, 5], "second": [6, 7, 8, 9, 10]}
        # y = {"first": [11, 12, 13, 14, 15], "second": [16, 17, np.nan, 19, 20]}
        # x = pd.DataFrame(x)
        # y = pd.DataFrame(y)
        # self.assertRaises(ValueError, lambda: bd_test(x=x, y=y))

        from sklearn.metrics.pairwise import euclidean_distances
        sigma = [[1, 0], [0, 1]]
        x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
        y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
        x = np.row_stack((x, y))
        dx = euclidean_distances(x, x)
        data_size = [50, 50]
        bd_value = bd_test(x=dx, size=data_size, dst=True)
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.059310879999)


if __name__ == '__main__':
    unittest.main()
    pass
