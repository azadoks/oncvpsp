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
 * running one TOML input, with the same sizes and no failure; every field
 * copied out, the same from both handles bit for bit, the test
 * configurations' energies the same as the text output's; the getters'
 * checks; the psp8, UPF and PSML files of both handles, the same byte for
 * byte; and a run of a broken input failing as an input failure at the
 * input stage. The stage handles: the all-electron atom alone, the same as
 * the run's; a run from it, every field the same as the whole run's; each
 * channel from it, the same as the run's convergence tables; and the model
 * core made again on the run's stages, without one the same as a run of the
 * input without one, and back to the input's the same as the first run.
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
    if (strcmp(group, "run") == 0) return strncmp(name, "channel.", 8) != 0;
    return strncmp(name, group, strlen(group)) == 0 && name[strlen(group)] == '.';
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

/* The stage handles of the input text, against the whole run hrun */
static void stages(void *hrun, int irel, const char *text, size_t len)
{
    oncv_dims_t dr, da;
    oncv_failure_t fl;
    void *ha = oncv_new(), *hh = oncv_new(), *hc = oncv_new();
    char *text0 = without_core(text);
    int ll, nn = 0, rank, shp[ONCV_MXRANK];

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
        check("channel", oncv_channel_toml(hc, ha, text, len, ll) == 0);
        same_slice(hc, "channel.cvgplt", hrun, "tests.convergence.cross", 3, ll, -1);
        same_slice(hc, "channel.ektail", hrun, "tests.convergence.tail", 2, ll, -1);
        same_slice(hc, "channel.optdg", hrun, "tests.convergence.opt", 2, ll, -1);
        same_slice(hc, "channel.wshare", hrun, "tests.convergence.wshare", 2, ll, -1);
        same_slice(hc, "channel.wcurve", hrun, "tests.weighted_cut.weight", 1, ll, -1);
        same_slice(hc, "channel.ecurve", hrun, "tests.convergence.curve", 3, ll, 1);
        oncv_get_int(hc, NAME("channel.l"), &nn, 1);
        check("channel.l", nn == ll);
    }
    check("channel from a run", oncv_channel_toml(hh, hrun, text, len, 0) == 0);
    oncv_channel_toml(hc, ha, text, len, 0);
    same_fields(hh, hc, "channel");
    ll = oncv_channel_toml(hc, ha, text, len, dr.lmax + 1);
    oncv_failure(hc, &fl);
    check("no such channel", ll % 100 == ONCV_FAIL_INPUT && fl.stage == ONCV_STAGE_INPUT);
    check("no ae in a failed handle", oncv_channel_toml(hh, hc, text, len, 0) % 100 ==
          ONCV_FAIL_INPUT);

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

    h1 = oncv_new();
    h2 = oncv_new();
    oncv_dims(h1, &d1);
    check("dims before a run", d1.version == ONCV_API_VERSION && d1.mmax == 0);
    oncv_set_output(h1, argv[3], strlen(argv[3]));
    ifail = oncv_run_toml(h1, irel, text, len, 0);
    check("run 1", ifail == 0);
    ifail = oncv_run_toml(h2, irel, text, len, 0);
    check("run 2", ifail == 0);
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

    ifail = oncv_run_toml(h2, irel, broken, strlen(broken), 0);
    oncv_failure(h2, &fl);
    check("broken input fails", ifail != 0);
    check("an input failure", fl.kind == ONCV_FAIL_INPUT && ifail % 100 == fl.kind);
    check("its message", fl.msg > 0 && ifail / 100 % 1000 == fl.msg);
    check("at the input", fl.stage == ONCV_STAGE_INPUT && fl.ll == -1);
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
