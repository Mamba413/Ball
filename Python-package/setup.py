# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:57:31 2018

@author: 99493
"""

from setuptools import setup, find_packages, Extension

cball_module = Extension('_cball',
                         sources=['src/cball_wrap.c', 'src/BD.c', 'src/utilities.c',
                                  'src/utilize_R.c', 'src/utilize_cross.c', 'src/cball.i'],
                         )

setup(name='Ball',
      version='0.1.0',
      author="Yue Hu",
      author_email="huyue45@mail2.sysu.edu.cn",
      maintainer="Jin Zhu",
      maintainer_email="zhuj37@mail2.sysu.edu.cn",
      description="Ball Python Package",
      long_description="""Hypothesis tests and sure independence screening (SIS) procedure based on ball statistics, 
      including ball divergence, ball covariance, and ball correlation, are developed to analyze complex data. 
      The ball divergence and ball covariance based distribution-free tests are implemented to examine equality of 
      multivariate distributions and independence between random vectors of arbitrary dimensions. Furthermore, 
      a generic non-parametric SIS procedure based on ball correlation and all of its variants are implemented to tackle 
      the challenge in the context of ultra high dimensional data.""",
      packages=find_packages(),
      install_requires=[
            'numpy',
            'scikit-learn'
      ],
      license="GPL-3",
      url="https://github.com/Mamba413/Ball",
      classifiers=[
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.6"
      ],
      ext_modules=[cball_module],
      )