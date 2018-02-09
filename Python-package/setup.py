# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:57:31 2018

@author: 99493
"""

from distutils.core import setup, Extension

cball_module = Extension('_cball',
                           sources=['cball_wrap.c', 'cball.c', 'BD.c', 'utilities.c'],
                           )
 
setup (name = 'cball',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [cball_module],
       py_modules = ["cball"],
       )