%module cball

%{
#define SWIG_FILE_WITH_INIT
#include "BD.h"
#include "BI.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

void bd_stat(double *, double *, int *, int *, int *, int *, int *, int *);
void bd_test(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *);
void bcov_stat(double *, double *, double *, int *, int *, int *, int *, int *);
void bcov_test(double *, double *, double *, double *, int *, int *, int *, int *, int *, int *);