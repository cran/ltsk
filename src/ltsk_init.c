#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
//extern void lk_main(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP RdistC(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"RdistC",  (DL_FUNC) &RdistC, 2},
  {NULL, NULL, 0}
};
/* static const R_CMethodDef CEntries[] = {
    {"lk_main", (DL_FUNC) &lk_main, 18},
    {NULL, NULL, 0}
};
 */
void R_init_ltsk(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
