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
 * The message itself is in the run's text output (oncv_set_output,
 * oncv_set_output_buffer).
 *
 * Beside a whole run (oncv_run_toml), a handle may hold a stage of one,
 * from which others start: the all-electron atom alone (oncv_ae_toml),
 * from which runs (oncv_run_ae_toml) and single channels
 * (oncv_channel_toml) start without solving it again; and a run's stages
 * before the model core, which every run keeps, from which oncv_coupled
 * makes the model core and the tests again with another model core. Each
 * call replaces what the handle held; the handle it starts from is read
 * only, and may be the handle itself.
 */

#ifndef ONCVPSP_H
#define ONCVPSP_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The version of this header; oncv_api_version() is the library's */
#define ONCV_API_VERSION 6

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
#define ONCV_STAGE_CONFINE 16

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

/* The library's build: its git commit, the state of the source against it
 * and the build time, a line each ("commit ", "source ", "compiled " and
 * the value), the first ntext bytes into text, not terminated, and its
 * whole length in *nout; ntext 0 gives the size to allocate */
int oncv_build_info(char *text, size_t ntext, size_t *nout);

/* A new run handle, its text output going to /dev/null */
void *oncv_new(void);

/* Frees a run handle and all it holds; NULL is ignored */
void oncv_free(void *handle);

/* The file the handle's runs write their text output to, npath bytes,
 * replaced at each run */
int oncv_set_output(void *handle, const char *path, size_t npath);

/* With on non-zero, the handle's runs write their text output to a buffer
 * the handle keeps, replaced at each run or stage, rather than to the file
 * of oncv_set_output; with on 0, to that file again */
int oncv_set_output_buffer(void *handle, int on);

/* The text output of the handle's last run or stage from its buffer: its
 * first ntext bytes into text, not terminated, and its whole length in
 * *nout (0 without a buffer); ntext 0 gives the size to allocate */
int oncv_output_text(void *handle, char *text, size_t ntext, size_t *nout);

/* The pseudopotential and its tests of the TOML input text (ntext bytes),
 * with the relativistic treatment irel, into the handle; skip_ld non-zero
 * leaves out the log derivatives */
int oncv_run_toml(void *handle, int irel, const char *text, size_t ntext,
                  int skip_ld);

/* The all-electron atom alone of the TOML input text, with the
 * relativistic treatment irel, into the handle: the fields ae.* */
int oncv_ae_toml(void *handle, int irel, const char *text, size_t ntext);

/* As oncv_run_toml, from the all-electron atom of ae_handle (oncv_ae_toml,
 * or any run's) rather than solving it; the input must be of that atom:
 * its irel, reference configuration, functional, mesh and speed of light
 * (else an input failure at the all-electron stage). The test
 * configurations may differ. The results are those of oncv_run_toml bit
 * for bit */
int oncv_run_ae_toml(void *handle, void *ae_handle, const char *text,
                     size_t ntext, int skip_ld);

/* The channel of angular momentum ll alone of the TOML input text, from
 * the all-electron atom of ae_handle (as oncv_run_ae_toml), into the
 * handle: the fields channel.*, the pseudo-waves, semi-local potential
 * and convergence curves as the whole run makes them. A first pseudo-wave
 * with a node is pspot's failure: ONCV_FAIL_CONSTRUCTION, msg 122
 * (msg_pspot_node). The mode is
 *   0  the channel;
 *   1  also the local potential (channel.vloc; for lloc < 4 from the
 *      local channel, optimized for it when ll is another) and the
 *      channel's projectors (channel.evkb, vkb, bbeig; none for the local
 *      channel), the same as the whole run's, whose failures are at
 *      ONCV_STAGE_VKB;
 *   2  the node test alone: the channel up to the node test of each
 *      kappa's first projector, which fails as above, as do the steps
 *      before it; the handle then holds no fields channel.*.
 * Another mode is ONCV_FAIL_INPUT, msg 185 (msg_run_channel_mode) */
int oncv_channel_toml(void *handle, void *ae_handle, const char *text,
                      size_t ntext, int ll, int mode);

/* The local potential alone of the TOML input text, from the all-electron
 * atom of ae_handle, into the handle: the fields local.* and, for
 * lloc < 4, its local channel (optimized with one projector), the same as
 * the whole run's; for oncv_channel_local_toml and oncv_assemble_toml */
int oncv_local_toml(void *handle, void *ae_handle, const char *text,
                    size_t ntext);

/* The channel ll as oncv_channel_toml, with the local potential of
 * local_handle (oncv_local_toml) in place of making it, if not NULL; it
 * must be of this input and atom (else ONCV_FAIL_INPUT, msg 187). With
 * tests non-zero (mode 1 only), also the channel's log derivatives, on
 * the window energies (emin, emax, de; Ha), NULL for the input's, and its
 * ghost tests (channel.logder.*, channel.ghosts.*). They are made on the
 * local potential before the ionic tail a whole run switches it to beyond
 * the largest core radius rcmax, which depends on every channel and the
 * model core: the log derivatives are the run's bit for bit when they read
 * the potential no further than rcmax (channel.logder.exact), near them
 * otherwise; the ghost tests reach beyond it, and are near the run's */
int oncv_channel_local_toml(void *handle, void *ae_handle,
                            void *local_handle, const char *text,
                            size_t ntext, int ll, int mode, int tests,
                            const double *energies);

/* The run of the TOML input text as oncv_run_ae_toml from the all-electron
 * atom of ae_handle, with the channels of the handles channels
 * (oncv_channel_toml or oncv_channel_local_toml, mode 0 or 1; nchannels
 * of them, in any order) in place of optimizing them: one of each l from
 * 0 to lmax, the local channel (lloc < 4) that of local_handle unless one
 * of the channels is of l = lloc. Each must be of this input and atom
 * (else ONCV_FAIL_INPUT, msg 186; a missing or doubled l msg 190). The
 * results are those of oncv_run_ae_toml bit for bit; the text output has
 * none of the channels' lines */
int oncv_assemble_toml(void *handle, void *ae_handle, void *local_handle,
                       void *const *channels, int nchannels, const char *text,
                       size_t ntext, int skip_ld);

/* The model core, the screening and the tests made again from the stages
 * of base's run (which a run keeps even if its model core failed), with
 * the model core icmod, fcfact and rcfact in place of the input's (fcfact
 * 0 for icmod 0), into the handle, which then holds those stages and the
 * results oncv_run_toml gives for the input with that model core */
int oncv_coupled(void *handle, void *base, int icmod, double fcfact,
                 double rcfact, int skip_ld);

/* The confined atom of the run of the handle run (oncv_run_toml,
 * oncv_run_ae_toml or oncv_coupled) into the handle, as oncvpsp_confine.x
 * makes it of the run's all-electron and pseudopotential data files: its
 * lines in the text output (oncv_set_output_buffer) and its summary, the
 * fields confine.*. The radii of the sphere (bohr), nradii of them, 0 for
 * the default (1.25, 1.5 and 2 times the largest core radius), less than
 * 0 for no sphere; the channels' l, nchannels of them, 0 for all four;
 * energies emin, emax and the scan's step (Ha), NULL for the input's
 * window of the log derivatives and 0.1; the walls' strengths (Ha/bohr**2),
 * nwalls of them, less than 0 for the default (0.01, 0.03 and 0.1), 0 for
 * none; r0 the wall's onset (bohr) and ecutmax the basis's largest cutoff
 * (Ha), 0 for the largest core radius and 200. The handle may be run */
int oncv_confine(void *handle, void *run, const double *radii, int nradii,
                 const int *channels, int nchannels, const double *energies,
                 const double *walls, int nwalls, double r0, double ecutmax);

/* Writes the pseudopotential of the handle's run as the file format
 * ("psp8", "upf" or "psml", nformat bytes) to path (npath bytes) */
int oncv_write_files(void *handle, const char *format, size_t nformat,
                     const char *path, size_t npath);

/* The failure of the handle's last run; kind 0 if it did not fail */
int oncv_failure(void *handle, oncv_failure_t *failure);

/* The sizes of the handle's results; without a finished run, irel, nkap,
 * iexc and mmax of its all-electron atom if it holds one, and all 0 but
 * version else */
int oncv_dims(void *handle, oncv_dims_t *dims);

/* The largest rank of a field */
#define ONCV_MXRANK 5

/* The results are copied out by field, named after the components of the
 * Fortran results (src/results_m.f90; pot.* src/pseudo_m.f90, tests.*
 * src/test_results_m.f90; core.* the model core's parameters,
 * src/results_m.f90; warnings.rec the warnings in the order written, 5
 * by nwarn (at most 64): the warning, its stage, l, ikap and index, as
 * failure_m's wrn_* lists them; ae.* the all-electron atom,
 * src/ae_state_m.f90, of an all-electron handle or a run's; local.* the
 * local potential of oncv_local_toml; confine.* the
 * confined atom's summary, src/confine_summary_m.f90, -1 for a part left
 * out, cut the cutoffs at 1e-4 and 1e-5 Ha of all the edges then of the
 * top edges, inf whether each is never reached; channel.* a
 * channel, src/channel_m.f90, its arrays per projector and kappa (mxprj,
 * nkap) and the convergence curves to the end of the longest; vloc, evkb,
 * vkb and bbeig of size 0 without mode 1), in Fortran
 * order: the first index fastest, so
 * a C caller indexes shp[0] innermost, and numpy reshapes with order='F'.
 * The sizes in the Fortran comments use oncv_dims and constants_m
 * (mxprj=5, mxl=3, mxstate=30, mxcnf=5, nrhomod=5, nerrlv=7, mxkap=2,
 * mxghost=20).
 * A field the run did not compute (esr without irel 2, the log
 * derivatives with skip_ld) has the size 0; a field of a stage the handle
 * does not hold (channel.* of a run, pot.* of an all-electron atom) is an
 * input failure. */

/* Field ifield (1, 2, ...): its kind, 'r' real or 'i' integer (a logical
 * as 0 or 1), a blank and its name, into name (at most nname bytes, not
 * terminated); returns the bytes written, 0 past the last field */
int oncv_field_name(int ifield, char *name, size_t nname);

/* What oncv_field_info gives of a field */
#define ONCV_INFO_AXES 1  /* its axes' names, blank-separated, first index
                           * first; an axis of different quantities as
                           * their names in brackets, "[ecut,err]" */
#define ONCV_INFO_UNITS 2 /* its units, empty for none; an axis of
                           * quantities' units in brackets, "[Ha,1]" */
#define ONCV_INFO_DOC 3   /* what it is, one sentence */

/* What field ifield is (what, ONCV_INFO_*): its first ntext bytes into
 * text, not terminated, and its whole length in *nout; ntext 0 gives the
 * size to allocate. The axes are r the radial mesh, iproj the projector,
 * l, ikap, istate the state, level the error level, iq the q grid, ie the
 * energy, icfg the test configuration, irow its row, ighost, deriv the
 * derivative, iv the valence state, iwarn the warning, icrad the basin
 * radius and ipow the power. An ifield past the last field or another
 * what is an input failure. */
int oncv_field_info(int ifield, int what, char *text, size_t ntext,
                    size_t *nout);

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
