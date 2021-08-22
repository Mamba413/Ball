# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:03:59 2019

@author: zhangyanhang
"""

import unittest
from unittest import TestCase
from Ball.bd import bd_test, bd
from Ball.bcorsis import bcorsis
from Ball.bcov import bcov_test, bcov
from Ball.wrap_c import bd_test_wrap_c, bcor_test_wrap_c, bcov_test_wrap_c, kbcov_test_wrap_c
import numpy as np
import math


class Test_bd(TestCase):
    def test_bd(self):
        np.random.seed(7654567)
        x = np.random.normal(0, 1, 50)
        y = np.random.normal(1, 1, 50)
        bd_value = bd(x, y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.196408479999)

        x = np.random.normal(0, 1, 100).reshape(50, 2)
        y = np.random.normal(3, 1, 100).reshape(50, 2)
        bd_value = bd(x, y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.5681075200000011)

        from sklearn.metrics.pairwise import euclidean_distances
        sigma = [[1, 0], [0, 1]]
        x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
        y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
        bd_value = bd(x, y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.13847583999999966)  # 0.059310879999

        n = 90
        x = np.random.normal(0, 1, n)
        bd_value = bd(x, size=np.array([40, 50]))
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.011700674999999966)
        
        x = [np.random.normal(0, 1, num) for num in [40, 50]]
        x = np.hstack(x)
        bd_value = bd(x, [40, 50])
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.9639094650205711)


    def test_bd_test(self):
        x = [1, 2, 3, 4, 5]
        y = [1, 2, 3, 4, 5]
        bd_value = bd_test(x, y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.0)

        np.random.seed(7654567)
        x = np.random.normal(0, 1, 50)
        y = np.random.normal(1, 1, 50)
        bd_value = bd_test(x, y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.196408479999)

        x = np.random.normal(0, 1, 100).reshape(50, 2)
        y = np.random.normal(3, 1, 100).reshape(50, 2)
        bd_value = bd_test(x, y)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.5681075200000011)


        x = np.random.normal(0, 1, 100).reshape(50, 2)
        y = np.random.normal(10, 1, 100).reshape(50, 2) 
        z = np.random.normal(100, 1, 100).reshape(50, 2) 
        bd_value = bd_test(x, y, z)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 2.0604000000000022)        
        
        bd_value = bd_test(x, y, z, weight = "max")
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 1.3736000000000015)         
        
        n = 90
        x = np.random.normal(0, 1, n)
        bd_value = bd_test(x, size=np.array([40, 50]))
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.009086599999999997)

        x = [np.random.normal(0, 1, num) for num in [40, 50]]
        x = np.hstack(x)
        bd_value = bd_test(x, [40, 50])
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.9639094650205713)

        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = np.array(x, dtype=np.double)
        bd_value = bd_test(x, size=np.array([5, 5]))
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.7231999999999997)

        x = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]
        bd_value = bd_test(x)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 2.403199999999999)

        from sklearn.metrics.pairwise import euclidean_distances
        sigma = [[1, 0], [0, 1]]
        x = np.random.multivariate_normal(mean=[0, 0], cov=sigma, size=50)
        y = np.random.multivariate_normal(mean=[1, 1], cov=sigma, size=50)
        x = np.row_stack((x, y))
        dx = euclidean_distances(x, x)
        data_size = [50, 50]
        bd_value = bd_test(dx, size=data_size, dst=True)
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.10779759999999977)
        
        
    def test_bcorsis(self):
        np.random.seed(1000)
        n = 150
        p = 3000
        mean = np.zeros(p)
        cov = np.array([0.5]*np.square(p)).reshape((p, p))
        cov[np.diag_indices(3000)] = 1
        x = np.random.multivariate_normal(mean, cov, n)
        error = np.random.normal(0, 1, n)
        y = 4*np.square(x[:, 2])+6*np.square(x[:, 1])+8*x[:, 3]-10*x[:,4]+error
        x_num = np.ones(3000)
        target = [4, 1, 924, 2, 692, 3, 400, 2241, 2839, 2194, 170]
        result = bcorsis(y, x, x_num, method="lm", params = [5, 3], d = 11)
        self.assertAlmostEqual(result, target)
        
        x = np.random.normal(0, 1, n*p).reshape((n, p))
        error = np.random.normal(0, 1, n)
        y = 3 * x[:, 1] * x[:, 5] * x[:, 10] + error
        target = [10, 5, 1, 1118, 555, 1174, 2361, 567, 1599, 1739]
        result = bcorsis(y, x, x_num, method = "interaction", d=10)
        self.assertAlmostEqual(result, target)
        
        y = 3 * x[:, 1] + 5 * np.square(x[:, 3]) + error
        target = [3, 1, 2607, 2374, 762]
        result = bcorsis(y, x, x_num, d=5, weight="prob")
        self.assertAlmostEqual(result, target)
        
        target = [3, 1, 762, 2374, 2607]
        result = bcorsis(y, x, x_num, d=5, weight="chisq")
        self.assertAlmostEqual(result, target)    
        
        
    def test_bcov_test(self):
        np.random.seed(2)
        error = np.random.uniform(-0.3, 0.3, 50)
        x = np.random.uniform(0, 4*math.pi, 50)
        y = np.cos(x) + error
        bcov_value = bcov_test(x, y)[0]
        self.assertAlmostEqual(bcov_value, 0.0021155347839999917)
        bcov_value = bcov_test(x, y, weight = "prob")[0]
        self.assertAlmostEqual(bcov_value, 0.05600363939468131)
        
        x = np.random.normal(0, 1, 50)
        y = [1 if i > 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
        z = [1 if i < 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
        bcov_value = bcov_test(x, y, z)[0]
        self.assertAlmostEqual(bcov_value, 0.0025298346891263934)
        
        x = np.random.normal(-math.pi, math.pi, 100).reshape((50, 2))
        error = np.random.uniform(-0.1, 0.1, 50)
        y = np.sin(x[:, 0] + x[:, 1]) + error
        bcov_value = bcov_test(x, y, weight = "prob")[0]
        self.assertAlmostEqual(bcov_value, 0.03201829726950482)
        bcov_value = bcov_test(x, y, weight = "chisq")[0]
        self.assertAlmostEqual(bcov_value, 0.05011815702428237)
        
    def test_bcov(self):
        np.random.seed(2)
        error = np.random.uniform(-0.3, 0.3, 50)
        x = np.random.uniform(0, 4*math.pi, 50)
        y = np.cos(x) + error
        bcov_value = bcov(x, y)[0]
        self.assertAlmostEqual(bcov_value, 0.0021155347839999917)
        bcov_value = bcov(x, y, weight = "prob")[0]
        self.assertAlmostEqual(bcov_value, 0.05600363939468131)
        
        x = np.random.normal(0, 1, 50)
        y = [1 if i > 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
        z = [1 if i < 0 else 0 for i in x] * x + np.random.normal(0, 1, 50)
        bcov_value = bcov(x, y, z)[0]
        self.assertAlmostEqual(bcov_value, 0.0025298346891263934)
        
        x = np.random.normal(-math.pi, math.pi, 100).reshape((50, 2))
        error = np.random.uniform(-0.1, 0.1, 50)
        y = np.sin(x[:, 0] + x[:, 1]) + error
        bcov_value = bcov(x, y, weight = "prob")[0]
        self.assertAlmostEqual(bcov_value, 0.03201829726950482)
        bcov_value = bcov(x, y, weight = "chisq")[0]
        self.assertAlmostEqual(bcov_value, 0.05011815702428237)
        
if __name__ == '__main__':
    unittest.main()
    pass