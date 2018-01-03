## Test environments
* local OS X install, R 3.4.3
* win-builder (realease)

## R CMD check results
* OS X enviorment:
0 errors | 0 warnings | 0 note


* Win-builder:
3 NOTEs
```
...
checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jin Zhu<zhuj37@mail2.sysu.edu.cn>'
...
** running examples for arch 'i386' ... [58s] NOTE
Examples with CPU or elapsed time > 10s
           user system elapsed
bcov.test 39.45   0.02   39.80
bcorsis   14.13   0.04   14.29
** running examples for arch 'x64' ... [70s] NOTE
Examples with CPU or elapsed time > 10s
           user system elapsed
bcov.test 53.45   0.02   53.41
bcorsis   12.43   0.06   12.62
```

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
