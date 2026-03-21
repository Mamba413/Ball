# -*- coding: utf-8 -*-

import os
import glob
import numpy as np
import pybind11
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

this_directory = os.path.abspath(os.path.dirname(__file__))

# Look for C sources in ../src/ (development) or ./src/ (sdist / CI wheel build)
root_src = os.path.abspath(os.path.join(this_directory, '..', 'src'))
if not os.path.isdir(root_src):
    root_src = os.path.abspath(os.path.join(this_directory, 'src'))

with open(os.path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Collect C sources from top-level src/
c_sources = glob.glob(os.path.join(root_src, '*.c'))


class BuildExt(build_ext):
    """Custom build_ext: adds -std=c++11 for .cpp on GCC/Clang; uses /O2 on MSVC."""

    def build_extensions(self):
        if self.compiler.compiler_type == 'msvc':
            for ext in self.extensions:
                ext.extra_compile_args = ['/O2']
        else:
            original_compile = self.compiler._compile

            def custom_compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
                postargs = list(extra_postargs or [])
                if src.endswith('.cpp'):
                    postargs.append('-std=c++11')
                return original_compile(obj, src, ext, cc_args, postargs, pp_opts)

            self.compiler._compile = custom_compile
        super().build_extensions()


ext_modules = [
    Extension(
        'Ball._cball',
        sources=['Ball/_bindings.cpp'] + c_sources,
        include_dirs=[root_src, np.get_include(), pybind11.get_include()],
        extra_compile_args=['-O3'],
    ),
]

setup(
    name='Ball',
    version='0.3.0',
    author="Jin Zhu, Xueqin Wang",
    author_email="zhuj37@mail2.sysu.edu.cn, wangxq88@mail.sysu.edu.cn",
    maintainer="Jin Zhu",
    maintainer_email="zhuj37@mail2.sysu.edu.cn",
    description="Ball: A Python Package for Detecting Distribution Difference and Association in Metric Spaces",
    long_description=long_description,
    long_description_content_type='text/x-rst',
    packages=find_packages(),
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    install_requires=[
        'numpy',
        'scikit-learn',
        'pygam',
    ],
    license="GPL-3",
    url="https://github.com/Mamba413/Ball",
    classifiers=[
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    python_requires='>=3.8',
)
