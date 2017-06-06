#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP themetagenomics_picrust_otu(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"themetagenomics_picrust_otu", (DL_FUNC) &themetagenomics_picrust_otu, 2},
  {NULL, NULL, 0}
};

void R_init_themetagenomics(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
