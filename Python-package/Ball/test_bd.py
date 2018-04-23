#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/3/8 16:36
# @Author  : Mamba
# @File    : test_bd.py


import unittest
from unittest import TestCase
from Ball.bd_test import bd_test, bd
import numpy as np


class Test_bd(TestCase):
    def test_bd(self):
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = np.array(x, dtype=np.double)
        bd_value = bd(x=x, size=np.array([5, 5]))
        bd_value = bd_value[0]
        self.assertAlmostEqual(bd_value, 0.7232)

    def test_bd_test(self):
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = np.array(x, dtype=np.double)
        bd_value = bd_test(x=x, size=np.array([5, 5]))
        bd_value = bd_value[0][0]
        self.assertAlmostEqual(bd_value, 0.7232)

        x = [1, 2, 3, None, 5, 6, 7, 8, None, 10]
        x = np.array(x, dtype=np.double)
        bd_value = bd_test(x=x, size=np.array([5, 5]))
        bd_value = bd_value[0][0]





if __name__ == '__main__':
    unittest.main()
    pass
