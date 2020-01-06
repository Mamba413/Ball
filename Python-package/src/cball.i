%module cball

%{
#include "BD.h"
#include "bcor.h"
#include "BI.h"
#include "kbcov.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

void bd_test(double *, double *, double *, int *, int *, int *, int *, int *, int *);
void bcor_test(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
void bcov_test(double *, double *, double *, double *, int *, int *, int *, int *);
void kbcov_test(double *, double *, double *, int *, int *, int *, int *, int *);