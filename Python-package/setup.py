# -*- coding: utf-8 -*-

import re
import os
import sys
import numpy as np
from setuptools import setup, find_packages, Extension
is_windows = sys.platform.startswith('win32')

os_type = 'MS_WIN64' if is_windows else 'LINUX'
split_str = '\\' if is_windows else '/'
path = sys.path
python_path = sys.prefix
version = sys.version.split('|')[0].strip().split('.')
del version[-1]
version = ''.join(version)

# if is_windows:
#       sys_command = 'configure.sh ' + python_path + ' ' + version;
#       os.system(sys_command)
#       pass

path1 = python_path + split_str + "include"
path2 = python_path + split_str + "libs"

default_compile_args = ["-DNDEBUG", "-fopenmp", "-O2", "-Wall", "-std=gnu99", "-mtune=generic"]
ext_libraries = []
if is_windows:
      default_compile_args.append("-D%s" % os_type)
      default_compile_args.append("-I%s" % path1)
      default_compile_args.append("-L%s" % path2)
      ext_libraries = ['vcruntime140']
      pass

cball_module = Extension('Ball._cball',
                         sources=['src/cball.i', 'src/BD.c', 'src/utilities.c', 'src/kbd.c', 'src/bcor.c', 'src/BI.c', 'src/kbcov.c'],
                         extra_compile_args= default_compile_args,
                         extra_link_args=['-lgomp'],
                         libraries=ext_libraries,
                         language='c',
                         include_dirs = ["Ball"]
                         )

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
      long_description = f.read()

setup(name='Ball',
      version='0.2.9',
      author="Jin Zhu, Xueqin Wang",
      author_email="zhuj37@mail2.sysu.edu.cn, wangxq88@mail.sysu.edu.cn",
      maintainer="Jin Zhu",
      maintainer_email="zhuj37@mail2.sysu.edu.cn",
      description="Ball: A Python Package for Detecting Distribution Difference and Association in Metric Spaces",
      long_description=long_description,
      long_description_content_type='text/x-rst',
      packages=find_packages(),
      install_requires=[
            'numpy',
            'scikit-learn',
            'pygam'
      ],
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
      ext_modules=[cball_module],
      )
