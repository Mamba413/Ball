Ball Statistics
================

Introdution
------------
The fundamental problems for data mining, statistical analysis, and machine learning are:
- whether several distributions are different?
- whether random variables are dependent?
- how to pick out useful variables/features from a high-dimensional data?

These issues can be tackled by using **bd_test**, **bcov_test**, and **bcorsis** functions in the **Ball** package, respectively. They enjoy following admirable advantages:
- available for most of datasets (e.g., traditional tabular data, brain shape, functional connectome, wind direction and so on)
- insensitive to outliers, distribution-free and model-free;
- theoretically guaranteed and computationally efficient.

Installation
------------
- Pypi version         

To install the **Ball** Python package from Pypi, just run:        
```
pip install Ball
```

- Building Ball library for Python for Windows with MinGW-w64 (Advanced)

You could download MinGW (https://sourceforge.net/projects/mingw/) and then
add the path ``MinGW/bin`` to system environment variable "path".
Anaconda3 is also in needed, and the version should be greater than 3.4. You should 
add all the related path of Anaconda3 to system environment variable "path",
as well as the path of ``MinGW/bin``.