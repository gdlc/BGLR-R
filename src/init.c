#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void C_read_bed(char **, int *, int *, int *, int *);
extern void C_write_bed(char **, int *, int *, int *);

extern SEXP C_read_ped(SEXP);
extern SEXP C_sample_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sample_beta_BB_BCp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sample_beta_BB_BCp_groups(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sample_beta_groups(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sample_beta_lower_tri(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static R_NativePrimitiveArgType C_read_bed_t[] = {
    STRSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType C_write_bed_t[] = {
    STRSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
    {"C_read_bed", (DL_FUNC) &C_read_bed, 5, C_read_bed_t},
    {"C_write_bed", (DL_FUNC) &C_write_bed, 4, C_write_bed_t},
    {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
    {"C_read_ped", (DL_FUNC) &C_read_ped, 1},
    {"C_sample_beta", (DL_FUNC) &C_sample_beta, 9},
    {"C_sample_beta_BB_BCp", (DL_FUNC) &C_sample_beta_BB_BCp, 11},
    {"C_sample_beta_BB_BCp_groups", (DL_FUNC) &C_sample_beta_BB_BCp_groups, 13},
    {"C_sample_beta_groups", (DL_FUNC) &C_sample_beta_groups, 11},
    {"C_sample_beta_lower_tri", (DL_FUNC) &C_sample_beta_lower_tri, 9},
    {NULL, NULL, 0}
};

void R_init_BGLR(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
