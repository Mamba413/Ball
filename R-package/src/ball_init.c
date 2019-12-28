#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bcor_test(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void bcov_test(double *, double *, double *, double *, int *, int *, int *, int *);
extern void kbcov_test(double *, double *, double *, int *, int *, int *, int *, int *);
extern void bd_test(double *, double *, double *, int *, int *, int *, int *, int *, int *);
extern void SRCT_new(double *, int *, int *, double *, int *, double *);
extern void bdd_matrix_bias(double *, double *, int *, int *);
extern void bdd_matrix_bias_two_group(double *, double *, int *, int *, int *);

static const R_CMethodDef CEntries[] = {
  {"bcor_test", (DL_FUNC) &bcor_test, 11},
  {"bcov_test", (DL_FUNC) &bcov_test, 8},
  {"kbcov_test", (DL_FUNC) &kbcov_test, 8},
  {"bd_test",   (DL_FUNC) &bd_test,   9},
  {"SRCT_new",  (DL_FUNC) &SRCT_new,  6},
  {"bdd_matrix_bias",           (DL_FUNC) &bdd_matrix_bias,           4},
  {"bdd_matrix_bias_two_group", (DL_FUNC) &bdd_matrix_bias_two_group, 5},
  {NULL, NULL, 0}
};

void R_init_Ball(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}