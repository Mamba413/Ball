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
- Pypi version (stable)         

To install the **Ball** Python package from Pypi, just run:        
```
pip install Ball
```

- Git version (development)
## Windows
Ball support compilation with MinGW. You need to install MinGW (https://sourceforge.net/projects/mingw/).
Next, you should clone or download the Ball repo. 
Then run the following from the root of the Ball directory: 
```
sh configure.sh your_python_path your_python_version_number # e.g.  sh configure.sh C:/anaconda3 35
```

Authorship
-----------
Jin Zhu (zhuj37@mail2.sysu.edu.cn), Xueqin Wang (wangxq88@mail2.sysu.edu.cn)