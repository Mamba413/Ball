#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/4/14 17:37
# @Author  : Mamba
# @Site    : 
# @File    : test_bcov.py

import unittest
from unittest import TestCase
from Ball.bcov import bcov_test, bcov, bcor
import numpy as np
import random
import math


class Test_bcov(TestCase):
    def test_bcov_test(self):
        random.seed(5)
        error = np.array([random.uniform(-0.3, 0.3) for _ in range(50)])
        x = np.array([random.uniform(0, 4 * math.pi) for _ in range(50)])
        y = np.cos(x) + error
        bcov_value = bcov_test(x, y)
        bcov_value = bcov_value[0][0]
        self.assertAlmostEqual(bcov_value, 0.00235927)

        random.seed(5)
        x = np.array([random.uniform(-math.pi, math.pi) for _ in range(100)]).reshape(50, 2)
        error = np.array([random.uniform(-0.3, 0.3) for _ in range(50)])
        y = np.sin([(x[i, 0]) ** 2 + x[i, 1] for i in range(0, 50)]) + error
        bcov_value = bcov_test(x=x, y=y)
        bcov_value = bcov_value[0][0]
        self.assertAlmostEqual(bcov_value, 0.00073367916)

        np.random.seed(5)
        x = np.random.normal(0, 1, 30)
        y = (x > 0) * x + np.random.normal(0, 1, 30)
        z = (x <= 0) * x + np.random.normal(0, 1, 30)
        data_list = [x.tolist(), y.tolist(), z.tolist()]
        kbcov_value = bcov_test(data_list)
        kbcov_value = kbcov_value[0][0]
        self.assertAlmostEqual(kbcov_value, 2.5328149917)

        x = np.array([1, 2, 3, 4, None, 6, 7, 8, 9, 10])
        self.assertRaises(ValueError, lambda: bcov_test(x))

        x = np.array([1, 2, 3, 4, None, 6, 7, 8, 9, 10])
        y = np.array([1, 3, 5])
        bcov_test(x, y)
        self.assertRaises(ValueError, lambda: bcov_test(x, y))

        x = np.array([1, 2, 3, 4, None, 6, 7, 8, 9, 10])
        y = np.array([5, 3, 2, 6, 5, 8, None, 1, 14, 7])
        self.assertRaises(ValueError, lambda: bcov_test(x, y))

    def test_bcov(self):
        np.random.seed(5)
        n = 50
        x = np.random.normal(0, 1, n)
        y = np.random.normal(0, 1, n)
        bcov_value = bcov(x, y)
        bcov_value = bcov_value[0]
        self.assertAlmostEqual(bcov_value, 0.00081812524)

    def test_bcor(self):
        np.random.seed(5)
        n = 50
        x = np.random.normal(0, 1, n)
        y = np.random.normal(0, 1, n)
        bcor_value = bcor(x, y)
        bcor_value = bcor_value[0]
        self.assertAlmostEqual(bcor_value, 0.0278907766)






if __name__ == '__main__':
    unittest.main()
    pass