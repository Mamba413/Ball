#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bcor_test(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void bcov_test(double *, double *, double *, double *, int *, int *, int *, int *);
extern void kbcov_test(double *, double *, double *, int *, int *, int *, int *, int *);
extern void bd_test(double *, double *, double *, int *, int *, int *, int *, int *, int *);
extern void SRCT_new(double *, int *, int *, double *, int *, double *);
extern void bd_gwas_screening(double *, double *, double *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void bd_gwas_refining_single(double *, double *, double *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void bdd_matrix_bias(double *, double *, int *, int *);
extern void bdd_matrix_bias_two_group(double *, double *, int *, int *, int *);
extern void sbd_C(double *, double *, int *, int *, int *, int *);

/* .Call calls */
extern SEXP _Ball_mq_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"bcor_test", (DL_FUNC) &bcor_test, 12},
  {"bcov_test", (DL_FUNC) &bcov_test, 8},
  {"kbcov_test", (DL_FUNC) &kbcov_test, 8},
  {"bd_test",   (DL_FUNC) &bd_test,   9},
  {"SRCT_new",  (DL_FUNC) &SRCT_new,  6},
  {"bd_gwas_refining_single", (DL_FUNC) &bd_gwas_refining_single, 14},
  {"bd_gwas_screening",       (DL_FUNC) &bd_gwas_screening,       14},
  {"bdd_matrix_bias",           (DL_FUNC) &bdd_matrix_bias,           4},
  {"bdd_matrix_bias_two_group", (DL_FUNC) &bdd_matrix_bias_two_group, 5},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_Ball_mq_cpp", (DL_FUNC) &_Ball_mq_cpp, 4},
    {NULL, NULL, 0}
};

void R_init_Ball(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}