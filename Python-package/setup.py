# -*- coding: utf-8 -*-
"""
"""

import re
import os
import sys
from setuptools import setup, find_packages, Extension
import numpy as np
os_type = 'MS_WIN64' if sys.platform.startswith('win32') else 'LINUX'
path = sys.path

pattern = re.compile(r'python')
length = []
for i in range(len(path)):
    temp = (path[i]).split("\\")
    if "Anaconda3" in temp or "Anaconda2" in temp or "anaconda3" in temp or "anaconda2" in temp:
        length.append(len(temp))
    else:
        length.append(100)
    if pattern.match(temp[-1]):
        version = "".join(re.findall(r"\d", temp[-1]))
ind = np.argmin(length)
python_path = path[ind]

path1 = "-I" + python_path + "\\include"
path2 = "-L" + python_path + "\\libs"

os.system('pre.sh ' + python_path + ' ' + version)

cball_module = Extension('Ball._cball',
                         sources=['src/cball_wrap.c', 'src/BD.c', 'src/utilities.c', 'src/kbd.c', 'src/bcor.c', 'src/BI.c', 'src/kbcov.c'],
                         extra_compile_args=["-DNDEBUG", "-fopenmp", "-O2", "-Wall", "-std=gnu99", "-mtune=generic", "-D%s" % os_type, path1, path2],
                         extra_link_args=['-lgomp'],
                         libraries=["vcruntime140"],
                         language='c',
                         include_dirs = ["Ball"]
                         )

setup(name='Ball',
      version='0.2.0',
      author="Xueqin Wang, Wenliang Pan, Heping Zhang, Hongtu Zhu, Yuan Tian, Weinan Xiao, Chengfeng Liu, Ruihuang Liu, Jin Zhu, Yanhang Zhang",
      author_email="zhangyh93@mail2.sysu.edu.cn",
      maintainer="Yanhang Zhang",
      maintainer_email="zhangyh93@mail2.sysu.edu.cn",
      description="Ball Python Package",
      long_description="""Hypothesis tests and sure independence screening (SIS) procedure based on ball statistics, including ball divergence , ball covariance , and ball correlation , 
                          are developed to analyze complex data in metric spaces, e.g, shape, directional, compositional and symmetric positive definite matrix data. The ball divergence 
                          and ball covariance based distribution-free tests are implemented to detecting distribution difference and association in metric spaces . Furthermore, several 
                          generic non-parametric feature selection procedures based on ball correlation, BCor-SIS and all of its variants, are implemented to tackle the challenge in the 
                          context of ultra high dimensional data.""",
      packages=find_packages(),
      install_requires=[
            'numpy',
            'scikit-learn',
            'pygam'
      ],
      license="GPL-3",
      url="https://github.com/Mamba413/Ball",
      classifiers=[
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7"
      ],
      python_requires='>=3.4',
      ext_modules=[cball_module],
      )
