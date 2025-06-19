# -*- coding: utf-8 -*-

import re
import os
import sys
import numpy as np
import subprocess
from setuptools import setup, find_packages

# Automatically build CFFI extension before installation
subprocess.check_call([sys.executable,
                       os.path.join(os.path.dirname(__file__), 'build_ffi.py')])


this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
      long_description = f.read()

setup(name='Ball',
      version='0.3.0',
      author="Jin Zhu, Xueqin Wang",
      author_email="zhuj37@mail2.sysu.edu.cn, wangxq88@mail.sysu.edu.cn",
      maintainer="Jin Zhu",
      maintainer_email="zhuj37@mail2.sysu.edu.cn",
      description="Ball: A Python Package for Detecting Distribution Difference and Association in Metric Spaces",
      long_description=long_description,
      long_description_content_type='text/x-rst',
      packages=find_packages(),
      setup_requires=['cffi'],
      install_requires=[
            'numpy',
            'scikit-learn',
            'cffi',
            'pygam'
      ],
      cffi_modules=['build_ffi.py:ffibuilder'],
      license="GPL-3",
      url="https://github.com/Mamba413/Ball",
      classifiers=[
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Information Analysis',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7'
      ],
      python_requires='>=3.5',
)
