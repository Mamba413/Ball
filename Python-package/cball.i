%module cball

%{
#define SWIG_FILE_WITH_INIT
#include "cball.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);

double bd_stat(double *, int, int, int, int);
void bd_test(double *, double *, double *, int, int, int, int, int, int);