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
 * Test of the C interface (include/oncvpsp.h): two handles open at once
 * running one TOML input, with the same sizes and no failure, the second's
 * text output into its buffer the same as the first's file; every field
 * copied out, the same from both handles bit for bit, the test
 * configurations' energies the same as the text output's; the getters'
 * checks; the psp8, UPF and PSML files of both handles, the same byte for
 * byte; and a run of a broken input failing as an input failure at the
 * input stage, its message in the buffer. The stage handles: the all-electron atom alone, the same as
 * the run's; a run from it, every field the same as the whole run's; each
 * channel from it, the same as the run's convergence tables, and its
 * projectors the same as the run's; inputs of another mesh or speed of
 * light refused (check_ae); and the model
 * core made again on the run's stages, without one the same as a run of the
 * input without one, and back to the input's the same as the first run.
 * The confined atom of the run: its summary with neither sphere nor walls,
 * and a radius inside the core radii, bad options and a handle without a
 * run refused.
 *
 * usage: test_c_api IREL INPUT.toml OUTPUT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "oncvpsp.h"

static int nfail = 0;

static void check(const char *label, int ok)
{
    if (!ok) {
        fprintf(stderr, "FAILED: %s\n", label);
        nfail++;
    }
}

#define NAME(s) (s), strlen(s)

/* whether a field is of the group: "ae" (ae.*), "channel" (channel.*) or
 * "run" (all but channel.*, which a run handle holds) */
static int in_group(const char *name, const char *group)
{
    if (strcmp(group, "run") == 0)
        return strncmp(name, "channel.", 8) != 0 && strncmp(name, "confine.", 8) != 0 &&
               strncmp(name, "local.", 6) != 0;
    return strncmp(name, group, strlen(group)) == 0 && name[strlen(group)] == '.';
}

/* The number of axes oncv_field_info gives field ifield, -1 without a
 * sentence on it */
static int info_rank(int ifield)
{
    char text[256];
    size_t nout = 0, kk;
    int rank = 0;
    if (oncv_field_info(ifield, ONCV_INFO_DOC, NULL, 0, &nout) != 0 || nout == 0)
        return -1;
    oncv_field_info(ifield, ONCV_INFO_AXES, NULL, 0, &nout);
    if (nout >= sizeof text) return -1;
    oncv_field_info(ifield, ONCV_INFO_AXES, text, nout, &nout);
    for (kk = 0; kk < nout; kk++)
        if (text[kk] != ' ' && (kk == 0 || text[kk - 1] == ' ')) rank++;
    return rank;
}

/* Every field of the group of the two handles: the shapes and values the
 * same */
static void same_fields(void *h1, void *h2, const char *group)
{
    char field[64];
    int ifield, nn, rank1, rank2, shp1[ONCV_MXRANK], shp2[ONCV_MXRANK], nseen = 0;
    size_t size, kk;
    for (ifield = 1; (nn = oncv_field_name(ifield, field, sizeof field - 1)) > 0;
         ifield++) {
        const char *name = field + 2;
        field[nn] = '\0';
        if (!in_group(name, group)) continue;
        nseen++;
        check("shape", oncv_shape(h1, NAME(name), &rank1, shp1) == 0);
        oncv_shape(h2, NAME(name), &rank2, shp2);
        check("shapes of two handles", rank1 == rank2 &&
              memcmp(shp1, shp2, sizeof(int) * (size_t)rank1) == 0);
        check(name, info_rank(ifield) == rank1);
        for (size = 1, kk = 0; kk < (size_t)rank1; kk++) size *= (size_t)shp1[kk];
        if (field[0] == 'r') {
            double *v1 = calloc(size + 1, sizeof *v1), *v2 = calloc(size + 1, sizeof *v2);
            check(name, oncv_get_double(h1, NAME(name), v1, size) == 0);
            oncv_get_double(h2, NAME(name), v2, size);
            check(name, memcmp(v1, v2, size * sizeof *v1) == 0);
            free(v1);
            free(v2);
        } else {
            int *v1 = calloc(size + 1, sizeof *v1), *v2 = calloc(size + 1, sizeof *v2);
            check(name, oncv_get_int(h1, NAME(name), v1, size) == 0);
            oncv_get_int(h2, NAME(name), v2, size);
            check(name, memcmp(v1, v2, size * sizeof *v1) == 0);
            free(v1);
            free(v2);
        }
    }
    check("fields", nseen > 10);
}

/* a real field into a new array, its shape in shp */
static double *get_double(void *handle, const char *name, int *rank, int *shp,
                          size_t *size)
{
    double *vv;
    int kk;
    if (oncv_shape(handle, NAME(name), rank, shp) != 0) return NULL;
    for (*size = 1, kk = 0; kk < *rank; kk++) *size *= (size_t)shp[kk];
    vv = calloc(*size + 1, sizeof *vv);
    if (oncv_get_double(handle, NAME(name), vv, *size) != 0) {
        free(vv);
        return NULL;
    }
    return vv;
}

/* a real field of two handles: the same shape and values */
static void same_field(void *h1, void *h2, const char *name)
{
    int rank1, rank2, shp1[ONCV_MXRANK], shp2[ONCV_MXRANK];
    size_t size1, size2;
    double *v1 = get_double(h1, name, &rank1, shp1, &size1);
    double *v2 = get_double(h2, name, &rank2, shp2, &size2);
    check(name, v1 && v2 && rank1 == rank2 && size1 == size2 &&
          memcmp(v1, v2, size1 * sizeof *v1) == 0);
    free(v1);
    free(v2);
}

/* a channel's field against the run's, whose index lpos is l: the same,
 * but along the channel's index nqpos (if >= 0), the curve's points, which
 * may end before the run's, whose points beyond are 0 */
static void same_slice(void *hc, const char *cname, void *hr, const char *rname,
                       int lpos, int ll, int nqpos)
{
    int crank, cshp[ONCV_MXRANK], rrank, rshp[ONCV_MXRANK], idx[ONCV_MXRANK];
    int kk, ok = 1;
    size_t csize, rsize, ir, ic, rest, stride;
    double *cv = get_double(hc, cname, &crank, cshp, &csize);
    double *rv = get_double(hr, rname, &rrank, rshp, &rsize);
    check(cname, cv != NULL && rv != NULL && rrank == crank + 1);
    if (!cv || !rv || rrank != crank + 1) goto done;
    for (kk = 0; kk < crank; kk++)
        if (kk != nqpos && cshp[kk] != rshp[kk < lpos ? kk : kk + 1]) ok = 0;
    check("the channel's shape as the run's", ok);
    for (ir = 0; ok && ir < rsize; ir++) {
        for (rest = ir, kk = 0; kk < rrank; kk++) {
            idx[kk] = (int)(rest % (size_t)rshp[kk]);
            rest /= (size_t)rshp[kk];
        }
        if (idx[lpos] != ll) continue;
        for (kk = lpos; kk < crank; kk++) idx[kk] = idx[kk + 1];
        if (nqpos >= 0 && idx[nqpos] >= cshp[nqpos]) {
            if (rv[ir] != 0.0) ok = 0;
            continue;
        }
        for (ic = 0, stride = 1, kk = 0; kk < crank; kk++) {
            ic += stride * (size_t)idx[kk];
            stride *= (size_t)cshp[kk];
        }
        if (cv[ic] != rv[ir]) ok = 0;
    }
    check(cname, ok);
done:
    free(cv);
    free(rv);
}

/* text with line inserted after the first at, or appended without at */
static char *with_line(const char *text, const char *at, const char *line)
{
    const char *pos = at ? strstr(text, at) : NULL;
    size_t nn = pos ? (size_t)(pos - text) + strlen(at) : strlen(text);
    char *out = malloc(strlen(text) + strlen(line) + 2);
    memcpy(out, text, nn);
    strcpy(out + nn, line);
    strcat(out, text + nn);
    return out;
}

/* the channel of the input other from the atom of ha refused with msg at
 * the atom's stage */
static void refused(void *hc, void *ha, char *other, const char *what, int msg)
{
    oncv_failure_t fl;
    int ifail = oncv_channel_toml(hc, ha, other, strlen(other), 0, 0);
    oncv_failure(hc, &fl);
    check(what, ifail % 100 == ONCV_FAIL_INPUT && fl.stage == ONCV_STAGE_AE &&
          fl.msg == msg);
    free(other);
}

/* the size of a field, -1 if it is none */
static long field_size(void *handle, const char *name)
{
    int rank, shp[ONCV_MXRANK], kk;
    long size = 1;
    if (oncv_shape(handle, NAME(name), &rank, shp) != 0) return -1;
    for (kk = 0; kk < rank; kk++) size *= shp[kk];
    return size;
}

/* text with the table [model_core_charge] replaced by icmod = 0 */
static char *without_core(const char *text)
{
    const char *head = "[model_core_charge]\n";
    const char *at = strstr(text, head), *end;
    char *out;
    size_t nn;
    if (!at) return NULL;
    end = strstr(at + strlen(head), "\n[");
    if (!end) end = text + strlen(text);
    else end++;
    nn = (size_t)(at - text);
    out = malloc(strlen(text) + 64);
    memcpy(out, text, nn);
    sprintf(out + nn, "%sicmod = 0\n\n%s", head, end);
    return out;
}

/* the value of the key of the table [model_core_charge], 0 without */
static double core_value(const char *text, const char *key)
{
    const char *at = strstr(text, "[model_core_charge]\n"), *end, *kk;
    double value = 0.0;
    if (!at) return 0.0;
    end = strstr(at + 1, "\n[");
    for (kk = strstr(at, key); kk && (!end || kk < end); kk = strstr(kk + 1, key))
        if (kk[-1] == '\n' && sscanf(kk + strlen(key), " = %lf", &value) == 1) break;
    return value;
}

/* The run assembled from channels made on the local potential alone, its
 * fields the whole run's bit for bit; each channel's log derivatives the
 * run's when exact */
static void assembled(void *hrun, void *ha, int irel, const char *text, size_t len)
{
    void *hl = oncv_new(), *hh = oncv_new(), *chans[ONCV_MXLIN + 1];
    oncv_dims_t dr;
    int ll, nn, nchan = 0, exact, rank, shp[ONCV_MXRANK];
    (void)irel;
    oncv_dims(hrun, &dr);
    check("local potential", oncv_local_toml(hl, ha, text, len) == 0);
    oncv_get_int(hl, NAME("local.lloc"), &nn, 1);
    check("local.lloc", nn == dr.lloc);
    for (ll = 0; ll <= dr.lmax; ll++) {
        void *hc = oncv_new();
        check("channel on the local potential",
              oncv_channel_local_toml(hc, ha, hl, text, len, ll, 1, 1, NULL) == 0);
        oncv_get_int(hc, NAME("channel.logder.exact"), &exact, 1);
        /* the kappas of the channel, 1 for l=0 of a fully relativistic run */
        same_slice(hc, "channel.logder.ae", hrun, "tests.logder.ae", 1, ll, 1);
        if (exact) same_slice(hc, "channel.logder.ps", hrun, "tests.logder.ps", 1, ll, 1);
        check("the tests need mode 1",
              oncv_channel_local_toml(hh, ha, hl, text, len, ll, 0, 1, NULL) % 100 ==
              ONCV_FAIL_INPUT);
        if (ll == dr.lloc) {
            oncv_free(hc);
            continue;
        }
        chans[nchan++] = hc;
    }
    check("assembled", oncv_assemble_toml(hh, ha, hl, chans, nchan, text, len, 0) == 0);
    same_fields(hh, hrun, "run");
    if (nchan > 0) {
        check("an l missing",
              oncv_assemble_toml(hh, ha, NULL, chans, nchan - (dr.lloc == 4), text, len, 0) %
              100 == ONCV_FAIL_INPUT);
    }
    check("no local potential", oncv_shape(hh, NAME("local.vloc"), &rank, shp) % 100 ==
          ONCV_FAIL_INPUT);
    while (nchan > 0) oncv_free(chans[--nchan]);
    oncv_free(hl);
    oncv_free(hh);
}

/* The stage handles of the input text, against the whole run hrun */
static void stages(void *hrun, int irel, const char *text, size_t len)
{
    oncv_dims_t dr, da;
    oncv_failure_t fl;
    void *ha = oncv_new(), *hh = oncv_new(), *hc = oncv_new();
    char *text0 = without_core(text);
    int ll, nn = 0, rank, shp[ONCV_MXRANK];
    double *vloc0 = NULL, *vloc;
    size_t nloc0 = 0, nloc;

    check("ae", oncv_ae_toml(ha, irel, text, len) == 0);
    oncv_dims(hrun, &dr);
    oncv_dims(ha, &da);
    check("ae dims", da.mmax == dr.mmax && da.irel == irel && da.lmax == 0);
    same_fields(ha, hrun, "ae");
    check("no channel of an ae handle", oncv_shape(ha, NAME("channel.pswf"), &rank, shp) %
          100 == ONCV_FAIL_INPUT);
    check("no results of an ae handle", oncv_shape(ha, NAME("pot.vkb"), &rank, shp) %
          100 == ONCV_FAIL_INPUT);

    check("run from ae", oncv_run_ae_toml(hh, ha, text, len, 0) == 0);
    same_fields(hh, hrun, "run");

    for (ll = 0; ll <= dr.lmax; ll++) {
        check("channel", oncv_channel_toml(hc, ha, text, len, ll, 0) == 0);
        check("no projectors in mode 0", field_size(hc, "channel.vkb") == 0 &&
              field_size(hc, "channel.vloc") == 0);
        oncv_channel_toml(hh, ha, text, len, ll, 1);
        same_field(hh, hc, "channel.pswf");
        same_field(hh, hc, "channel.ecurve");
        if (ll == dr.lloc) {
            check("no projectors of the local channel",
                  field_size(hh, "channel.vkb") == 0 && field_size(hh, "channel.evkb") == 0);
        } else {
            same_slice(hh, "channel.evkb", hrun, "pot.evkb", 1, ll, -1);
            same_slice(hh, "channel.vkb", hrun, "pot.vkb", 2, ll, -1);
            same_slice(hh, "channel.bbeig", hrun, "tests.bbeig", 1, ll, -1);
        }
        vloc = get_double(hh, "channel.vloc", &rank, shp, &nloc);
        check("channel.vloc", vloc != NULL && nloc == (size_t)dr.mmax);
        if (!vloc0) {
            vloc0 = vloc;
            nloc0 = nloc;
        } else {
            check("the same local potential from each channel", vloc && nloc == nloc0 &&
                  memcmp(vloc, vloc0, nloc * sizeof *vloc) == 0);
            free(vloc);
        }
        same_slice(hc, "channel.cvgplt", hrun, "tests.convergence.cross", 3, ll, -1);
        same_slice(hc, "channel.ektail", hrun, "tests.convergence.tail", 2, ll, -1);
        same_slice(hc, "channel.optdg", hrun, "tests.convergence.opt", 2, ll, -1);
        same_slice(hc, "channel.wshare", hrun, "tests.convergence.wshare", 2, ll, -1);
        same_slice(hc, "channel.wcurve", hrun, "tests.weighted_cut.weight", 1, ll, -1);
        same_slice(hc, "channel.ecurve", hrun, "tests.convergence.curve", 3, ll, 1);
        oncv_get_int(hc, NAME("channel.l"), &nn, 1);
        check("channel.l", nn == ll);
    }
    free(vloc0);
    assembled(hrun, ha, irel, text, len);
    for (ll = 0; ll <= dr.lmax; ll++) {
        check("the node test", oncv_channel_toml(hh, ha, text, len, ll, 2) == 0);
        check("no channel after the node test",
              oncv_shape(hh, NAME("channel.pswf"), &rank, shp) % 100 == ONCV_FAIL_INPUT);
    }
    ll = oncv_channel_toml(hh, ha, text, len, 0, 3);
    oncv_failure(hh, &fl);
    check("no mode 3", ll % 100 == ONCV_FAIL_INPUT && fl.msg == 185);
    check("channel from a run", oncv_channel_toml(hh, hrun, text, len, 0, 1) == 0);
    oncv_channel_toml(hc, ha, text, len, 0, 1);
    same_fields(hh, hc, "channel");
    ll = oncv_channel_toml(hc, ha, text, len, dr.lmax + 1, 0);
    oncv_failure(hc, &fl);
    check("no such channel", ll % 100 == ONCV_FAIL_INPUT && fl.stage == ONCV_STAGE_INPUT);
    check("no ae in a failed handle", oncv_channel_toml(hh, hc, text, len, 0, 0) % 100 ==
          ONCV_FAIL_INPUT);
    refused(hc, ha, with_line(text, NULL, "\n[log_mesh]\namesh = 1.02\n"),
            "another mesh", 177);
    refused(hc, ha, with_line(text, "[oncvpsp]\n", "c_light = 100.0\n"),
            "another speed of light", 178);

    check("toml without the model core", text0 != NULL);
    if (text0) {
        check("run without the model core",
              oncv_run_toml(hc, irel, text0, strlen(text0), 0) == 0);
        check("coupled without the model core",
              oncv_coupled(hh, hrun, 0, 0.0, 0.0, 0) == 0);
        same_fields(hh, hc, "run");
        check("coupled again on its own stages",
              oncv_coupled(hh, hh, (int)core_value(text, "icmod"),
                           core_value(text, "fcfact"), core_value(text, "rcfact"), 0) == 0);
        same_fields(hh, hrun, "run");
    }
    check("no stages in an ae handle", oncv_coupled(hh, ha, 0, 0.0, 0.0, 0) % 100 ==
          ONCV_FAIL_INPUT);
    check("run from its own ae", oncv_run_ae_toml(ha, ha, text, len, 0) == 0);
    same_fields(ha, hrun, "run");

    free(text0);
    oncv_free(ha);
    oncv_free(hh);
    oncv_free(hc);
}

/* the handle's buffered text output, terminated, in a new array */
static char *output_text(void *handle)
{
    size_t nout = 0, nagain = 0;
    char *text;
    oncv_output_text(handle, NULL, 0, &nout);
    text = malloc(nout + 1);
    oncv_output_text(handle, text, nout, &nagain);
    check("the text's size", nagain == nout);
    text[nout] = '\0';
    return text;
}

/* The confined atom of the run hrun, its parts left out, and its refusals
 * (oncvpsp_confine.x's own checks are scripts/c_api_confine.py's) */
static void confine(void *hrun)
{
    oncv_failure_t fl;
    void *hc = oncv_new();
    int nn = 0, ll = 7, rank, shp[ONCV_MXRANK];
    double rr = 0.1, de = -1.0, ee[3] = {1.0, -1.0, 0.1};

    check("confine without sphere or walls",
          oncv_confine(hc, hrun, NULL, -1, NULL, 0, NULL, NULL, 0, 0.0, 0.0) == 0);
    oncv_get_int(hc, NAME("confine.nmiss"), &nn, 1);
    oncv_get_double(hc, NAME("confine.max_de"), &de, 1);
    check("the parts left out", nn == -1 && de == -1.0);
    check("no run fields of a confine handle", oncv_shape(hc, NAME("pot.rr"), &rank,
          shp) % 100 == ONCV_FAIL_INPUT);
    ll = oncv_confine(hc, hrun, &rr, 1, NULL, 0, NULL, NULL, 0, 0.0, 0.0);
    oncv_failure(hc, &fl);
    check("a radius inside the core radii", ll % 100 == ONCV_FAIL_INPUT &&
          fl.stage == ONCV_STAGE_CONFINE && fl.msg == 182);
    ll = 7;
    check("a channel past 3", oncv_confine(hc, hrun, NULL, -1, &ll, 1, NULL, NULL, 0,
          0.0, 0.0) / 100 % 1000 == 184);
    check("a window upside down", oncv_confine(hc, hrun, NULL, -1, NULL, 0, ee, NULL, 0,
          0.0, 0.0) / 100 % 1000 == 184);
    check("no run", oncv_confine(hc, hc, NULL, -1, NULL, 0, NULL, NULL, 0, 0.0, 0.0) /
          100 % 1000 == 183);
    oncv_free(hc);
}

/* The test configurations' energies against the text output's test_config
 * lines (diagout, e24.16: 17 digits, which give the double back) */
static void same_as_text(void *handle, const char *path)
{
    double de[3 * 5], aa[3];
    char line[512];
    int ncfg, icfg, nseen = 0;
    FILE *fp = fopen(path, "r");
    check("text output", fp != NULL);
    if (!fp) return;
    oncv_get_int(handle, NAME("tests.configs.ncfg"), &ncfg, 1);
    check("configs.de", oncv_get_double(handle, NAME("tests.configs.de"), de, 15) == 0);
    while (fgets(line, sizeof line, fp)) {
        if (sscanf(line, "test_config %d %lf %lf %lf", &icfg, aa, aa + 1, aa + 2) != 4)
            continue;
        check("config in range", icfg >= 0 && icfg < 5);
        check("test_config as the text", aa[0] == de[3 * icfg] &&
              aa[1] == de[3 * icfg + 1] && aa[2] == de[3 * icfg + 2]);
        nseen++;
    }
    fclose(fp);
    check("test_config lines", nseen == ncfg && ncfg > 0);
}

static char *read_file(const char *path, size_t *len);

/* the files of the two handles' runs (OUTPUT.psp8.1, ...), the same byte
 * for byte, each beginning as its format does */
static void same_files(void *h1, void *h2, const char *out)
{
    const char *formats[3] = {"psp8", "upf", "psml"};
    const char *starts[3] = {"", "<UPF version", "<?xml"};
    char path1[4096], path2[4096];
    char *f1, *f2;
    size_t n1, n2;
    int ii;
    for (ii = 0; ii < 3; ii++) {
        snprintf(path1, sizeof path1, "%s.%s.1", out, formats[ii]);
        snprintf(path2, sizeof path2, "%s.%s.2", out, formats[ii]);
        check(formats[ii], oncv_write_files(h1, formats[ii], strlen(formats[ii]),
                                            path1, strlen(path1)) == 0 &&
                           oncv_write_files(h2, formats[ii], strlen(formats[ii]),
                                            path2, strlen(path2)) == 0);
        f1 = read_file(path1, &n1);
        f2 = read_file(path2, &n2);
        check("the files of two handles",
              f1 && f2 && n1 > 0 && n1 == n2 && memcmp(f1, f2, n1) == 0);
        check("the file's start", f1 && strncmp(f1, starts[ii], strlen(starts[ii])) == 0);
        free(f1);
        free(f2);
    }
    check("no such format", oncv_write_files(h1, "cpi", 3, path1, strlen(path1)) % 100 ==
          ONCV_FAIL_INPUT);
}

static char *read_file(const char *path, size_t *len)
{
    FILE *fp = fopen(path, "rb");
    char *text;
    long size;
    if (!fp) return NULL;
    fseek(fp, 0, SEEK_END);
    size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    text = malloc((size_t)size + 1);
    if (text && fread(text, 1, (size_t)size, fp) != (size_t)size) {
        free(text);
        text = NULL;
    }
    fclose(fp);
    if (text) {
        text[size] = '\0';
        *len = (size_t)size;
    }
    return text;
}

int main(int argc, char **argv)
{
    oncv_dims_t d1, d2;
    oncv_failure_t fl;
    const char *broken = "[oncvpsp]\natsym = 'O'\n";
    void *h1, *h2;
    size_t len;
    char *text;
    int irel, ifail, l1;

    if (argc != 4) {
        fprintf(stderr, "usage: test_c_api IREL INPUT.toml OUTPUT\n");
        return 2;
    }
    irel = atoi(argv[1]);
    text = read_file(argv[2], &len);
    if (!text) {
        fprintf(stderr, "test_c_api: cannot read %s\n", argv[2]);
        return 2;
    }
    check("api version", oncv_api_version() == ONCV_API_VERSION);
    {
        size_t nout = 0;
        char info[512];
        oncv_build_info(NULL, 0, &nout);
        check("build info", nout > 0 && nout < sizeof info);
        oncv_build_info(info, nout, &nout);
        check("its lines", strncmp(info, "commit ", 7) == 0 && memchr(info, '\n', nout));
    }

    {
        size_t nout = 0;
        check("field info past the last field",
              oncv_field_info(100000, ONCV_INFO_DOC, NULL, 0, &nout) != 0 && nout == 0);
        check("field info of no such kind", oncv_field_info(1, 4, NULL, 0, &nout) != 0);
    }

    h1 = oncv_new();
    h2 = oncv_new();
    oncv_dims(h1, &d1);
    check("dims before a run", d1.version == ONCV_API_VERSION && d1.mmax == 0);
    oncv_set_output(h1, argv[3], strlen(argv[3]));
    ifail = oncv_run_toml(h1, irel, text, len, 0);
    check("run 1", ifail == 0);
    oncv_set_output_buffer(h2, 1);
    ifail = oncv_run_toml(h2, irel, text, len, 0);
    check("run 2", ifail == 0);
    {
        size_t nfile;
        char *file = read_file(argv[3], &nfile), *buffer = output_text(h2);
        check("the buffer as the file", file && nfile > 0 && strlen(buffer) == nfile &&
              memcmp(file, buffer, nfile) == 0);
        free(file);
        free(buffer);
    }
    oncv_failure(h1, &fl);
    check("no failure", fl.kind == ONCV_FAIL_NONE && fl.msg == 0);
    oncv_dims(h1, &d1);
    oncv_dims(h2, &d2);
    check("dims", d1.mmax > 0 && d1.irel == irel && d1.lmax >= 0);
    check("dims of two handles", memcmp(&d1, &d2, sizeof d1) == 0);
    for (l1 = 0; l1 <= d1.lmax; l1++)
        check("projectors", d1.lloc == l1 || d1.nproj[l1] > 0);
    same_fields(h1, h2, "run");
    same_as_text(h1, argv[3]);
    same_files(h1, h2, argv[3]);
    {
        int rank, shp[ONCV_MXRANK], ival;
        double xx[2];
        check("rr's shape", oncv_shape(h1, NAME("pot.rr"), &rank, shp) == 0 &&
              rank == 1 && shp[0] == d1.mmax);
        check("a scalar", oncv_shape(h1, NAME("pot.zion"), &rank, shp) == 0 && rank == 0);
        check("no such field", oncv_shape(h1, NAME("pot.nothing"), &rank, shp) % 100 ==
              ONCV_FAIL_INPUT);
        check("wrong size", oncv_get_double(h1, NAME("pot.zion"), xx, 2) % 100 ==
              ONCV_FAIL_INPUT);
        check("wrong kind", oncv_get_int(h1, NAME("pot.zion"), &ival, 1) % 100 ==
              ONCV_FAIL_INPUT);
    }

    stages(h1, irel, text, len);
    confine(h1);

    ifail = oncv_run_toml(h2, irel, broken, strlen(broken), 0);
    oncv_failure(h2, &fl);
    check("broken input fails", ifail != 0);
    check("an input failure", fl.kind == ONCV_FAIL_INPUT && ifail % 100 == fl.kind);
    check("its message", fl.msg > 0 && ifail / 100 % 1000 == fl.msg);
    check("at the input", fl.stage == ONCV_STAGE_INPUT && fl.ll == -1);
    {
        char *buffer = output_text(h2);
        check("its message in the buffer", strstr(buffer, "ERROR") != NULL);
        free(buffer);
    }
    oncv_set_output_buffer(h2, 0);
    oncv_run_toml(h2, irel, broken, strlen(broken), 0);
    {
        size_t nout = 1;
        oncv_output_text(h2, NULL, 0, &nout);
        check("no buffer after it is off", nout == 0);
    }
    oncv_dims(h2, &d2);
    check("no results after a failure", d2.mmax == 0);
    oncv_dims(h1, &d1);
    check("the other handle keeps its results", d1.mmax > 0);

    oncv_free(h1);
    oncv_free(h2);
    oncv_free(NULL);
    free(text);
    if (nfail) {
        fprintf(stderr, "c_api: %d check(s) failed\n", nfail);
        return 1;
    }
    fprintf(stderr, "c_api: all checks passed\n");
    return 0;
}
