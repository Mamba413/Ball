## Resubmission
This is a resubmission. In this version I have:

* remove time consuming examples

## Test environments
* local OS X install, R 3.4.3
* Ubuntu 14.04 (on travis-ci), R 3.4.2
* win-builder (realease)

## R CMD check results
* OS X enviorment:
0 errors | 0 warnings | 0 note

* Ubuntu 14.04 enviorment:
0 errors | 0 warnings | 0 note

* Win-builder:
There was 1 NOTE:

```
checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jin Zhu <zhuj37@mail2.sysu.edu.cn>'

New submission
```

I think [this NOTE](https://stackoverflow.com/questions/36701433/r-package-building-with-devtoolsbuild-win-version-contains-large-components) is ignorable. It occurs because it is my first submission of the package.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
