/*  File src/init.c in package ergm.userterms, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2012-2019 Statnet Commons
 */
#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*The package's .Call entry points.  Declared here rather than in a header
  because the registration table below is their only other reference, and it
  is the one place their arity is recorded.*/
extern SEXP simEGP_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                     SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EGPHazard_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

/*Registering the .Call routines makes R check the argument count at call
  time.  Without it a mismatch between the R and C sides -- an argument added
  to one and not the other, or a stale shared library left over from an
  earlier build -- is undefined behaviour rather than an error, which is a
  singularly unpleasant thing to debug.  Keep the counts here in step with the
  signatures above.*/
static const R_CallMethodDef CallEntries[] = {
    {"simEGP_R",    (DL_FUNC) &simEGP_R,    18},
    {"EGPHazard_R", (DL_FUNC) &EGPHazard_R,  8},
    {NULL, NULL, 0}
};

/*This must be named for the package, or R never calls it: it was previously
  R_init_ergm_ego, inherited from the template this file was copied from, so
  the registration below had no effect at all.*/
void R_init_ergmgp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    /*Left TRUE deliberately: .Call() is used by name throughout the R code.*/
    R_useDynamicSymbols(dll, TRUE);
}
