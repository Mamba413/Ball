## Maintainer comments    

- We have removed the warning message when compiling C code on MacOS.

- The DOI and URL in the DESCRIPTION/README.md correspond to an AOS publication. 
Actually, the DOI and URL are definitely right. People can search the DOI through  search.crossref.org, and visit the website according to the URL. 

## Test environments
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Win-builder (devel)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Local computer: Darwin17.0, x86_64

## R CMD check results
NOTE (only in Win-builder, devel): 

```
Found the following (possibly) invalid URLs:
  URL: https://projecteuclid.org/euclid.aos/1525313077
    From: README.md
    Status: 500
    Message: Internal Server Error
    
Found the following (possibly) invalid DOIs:
  DOI: 10.1214/17-AOS1579
    From: DESCRIPTION
    Status: Internal Server Error
    Message: 500
```

## Reverse dependencies
There are no reverse dependencies.