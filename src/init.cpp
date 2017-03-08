#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "include/RMBA.h"

static R_CallMethodDef CallEntries[] = {
    {"MBAPoints", (DL_FUNC) &MBAPoints, 7},
    {"MBASurf", (DL_FUNC) &MBASurf, 12},
    {NULL, NULL, 0}
};

void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_sp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

