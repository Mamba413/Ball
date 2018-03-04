#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bcov_stat(double *, double *, double *, int *, int *, int *, int *);
extern void bcov_test(double *, double *, double *, double *, int *, int *, int *, int *, int *);
extern void bd_stat(double *, double *, int *, int *, int *, int *, int *);
extern void bd_test(double *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void kbd_stat(double *, double *, double *, int *, int *, int *, int *, int *);
extern void kbd_test(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *);
extern void SRCT(double *, double *, double *, double *, int *, double *);

static const R_CMethodDef CEntries[] = {
  {"bcov_stat", (DL_FUNC) &bcov_stat, 7},
  {"bcov_test", (DL_FUNC) &bcov_test, 9},
  {"bd_stat",   (DL_FUNC) &bd_stat,   7},
  {"bd_test",   (DL_FUNC) &bd_test,   10},
  {"kbd_stat",  (DL_FUNC) &kbd_stat,  8},
  {"kbd_test",  (DL_FUNC) &kbd_test,  10},
  {"SRCT",      (DL_FUNC) &SRCT,      6},
  {NULL, NULL, 0}
};

void R_init_Ball(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}