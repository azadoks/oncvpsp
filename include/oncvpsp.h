/*
 * Copyright (c) 2026 by Austin Zadoks
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * The C interface of ONCVPSP (liboncvpsp, src/oncvpsp_c_m.F90).
 *
 * A run handle holds one calculation: its input, given as TOML text, its
 * results and its failure. Handles share nothing, so two handles are two
 * independent calculations, in one thread or in several; one handle must
 * not be used by two threads at once. Every handle passed must come from
 * oncv_new and not yet be freed.
 *
 * The functions that return int return 0, or the failure's ifail: the
 * kind of failure in its last two digits (ONCV_FAIL_*), which message
 * above them, and the solver status above that; oncv_failure splits it,
 * with where the run was. The messages (msg) are the msg_* of
 * src/failure_m.f90, whose values are kept from one version to the next.
 * The message itself is in the run's text output (oncv_set_output).
 */

#ifndef ONCVPSP_H
#define ONCVPSP_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The version of this header; oncv_api_version() is the library's */
#define ONCV_API_VERSION 1

/* The kinds of failure (failure_m's fail_*) */
#define ONCV_FAIL_NONE 0
#define ONCV_FAIL_INPUT 1
#define ONCV_FAIL_CONVERGENCE 2
#define ONCV_FAIL_NOT_FOUND 3
#define ONCV_FAIL_CONSTRUCTION 4
#define ONCV_FAIL_NUMERICAL 5
#define ONCV_FAIL_IO 6
#define ONCV_FAIL_INTERNAL 7

/* The stages of a run (failure_m's stage_*) */
#define ONCV_STAGE_NONE 0
#define ONCV_STAGE_INPUT 1
#define ONCV_STAGE_AE 2
#define ONCV_STAGE_CHANNEL 3
#define ONCV_STAGE_VKB 4
#define ONCV_STAGE_VALENCE 5
#define ONCV_STAGE_MODEL_CORE 6
#define ONCV_STAGE_SCREENING 7
#define ONCV_STAGE_DIAG 8
#define ONCV_STAGE_GHOSTS 9
#define ONCV_STAGE_UNSCREEN 10
#define ONCV_STAGE_CONFIGS 11
#define ONCV_STAGE_PLOT 12
#define ONCV_STAGE_PHSFT 13
#define ONCV_STAGE_TESTS 14
#define ONCV_STAGE_OUTPUT 15

/* The largest l of the input's channels, l = 0..ONCV_MXLIN */
#define ONCV_MXLIN 5

/* A failure (failure_m's failure_t): 0 where a field does not apply, ll -1 */
typedef struct {
    int kind;  /* the kind of failure, ONCV_FAIL_* */
    int msg;   /* which message, a msg_* of failure_m, 0 if none */
    int ierr;  /* the solver status the failure stopped on, 0 if none */
    int stage; /* the stage of the run, ONCV_STAGE_* */
    int ll;    /* the channel's angular momentum, -1 if none */
    int ikap;  /* the spin-orbit index (1 or 2), 0 if none */
    int iproj; /* the projector, 0 if none */
} oncv_failure_t;

/* The sizes of a run's results (oncvpsp_c_m's dims_t) */
typedef struct {
    int version;               /* the library's ONCV_API_VERSION */
    int irel;                  /* 0 non-, 1 scalar-, 2 fully-relativistic */
    int nkap;                  /* kappas per l: 1, or 2 for irel 2 */
    int iexc;                  /* the functional */
    int mmax;                  /* size of the log radial mesh */
    int lmax;                  /* the largest l of the projectors */
    int lloc;                  /* the local potential's l, 4: polynomial */
    int nproj[ONCV_MXLIN + 1]; /* the projectors of each l */
    int nrl;                   /* points of the linear output mesh */
} oncv_dims_t;

/* The version of the library's interface */
int oncv_api_version(void);

/* A new run handle, its text output going to /dev/null */
void *oncv_new(void);

/* Frees a run handle and all it holds; NULL is ignored */
void oncv_free(void *handle);

/* The file the handle's runs write their text output to, npath bytes,
 * replaced at each run */
int oncv_set_output(void *handle, const char *path, size_t npath);

/* The pseudopotential and its tests of the TOML input text (ntext bytes),
 * with the relativistic treatment irel, into the handle; skip_ld non-zero
 * leaves out the log derivatives */
int oncv_run_toml(void *handle, int irel, const char *text, size_t ntext,
                  int skip_ld);

/* The failure of the handle's last run; kind 0 if it did not fail */
int oncv_failure(void *handle, oncv_failure_t *failure);

/* The sizes of the handle's results; all 0 but version without a
 * finished run */
int oncv_dims(void *handle, oncv_dims_t *dims);

/* The largest rank of a field */
#define ONCV_MXRANK 5

/* The results are copied out by field, named after the components of the
 * Fortran results (src/results_m.f90; pot.* src/pseudo_m.f90, tests.*
 * src/test_results_m.f90), in Fortran order: the first index fastest, so
 * a C caller indexes shp[0] innermost, and numpy reshapes with order='F'.
 * The sizes in the Fortran comments use oncv_dims and constants_m
 * (mxprj=5, mxl=3, mxstate=30, mxcnf=5, nrhomod=5, nerrlv=7, mxkap=2).
 * A field the run did not compute (esr without irel 2, the log
 * derivatives with skip_ld) has the size 0. */

/* Field ifield (1, 2, ...): its kind, 'r' real or 'i' integer (a logical
 * as 0 or 1), a blank and its name, into name (at most nname bytes, not
 * terminated); returns the bytes written, 0 past the last field */
int oncv_field_name(int ifield, char *name, size_t nname);

/* The rank of a field and its shape, shp[0..rank-1] */
int oncv_shape(void *handle, const char *name, size_t nname, int *rank,
               int shp[ONCV_MXRANK]);

/* The values of a real field, nvalues of them (its size, oncv_shape) */
int oncv_get_double(void *handle, const char *name, size_t nname,
                    double *values, size_t nvalues);

/* The values of an integer field, nvalues of them */
int oncv_get_int(void *handle, const char *name, size_t nname, int *values,
                 size_t nvalues);

#ifdef __cplusplus
}
#endif

#endif
