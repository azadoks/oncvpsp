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
 * input stage.
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

/* Every field of the two handles: the shapes and values the same */
static void same_fields(void *h1, void *h2)
{
    char field[64];
    int ifield, nn, rank1, rank2, shp1[ONCV_MXRANK], shp2[ONCV_MXRANK];
    size_t size, kk;
    for (ifield = 1; (nn = oncv_field_name(ifield, field, sizeof field - 1)) > 0;
         ifield++) {
        const char *name = field + 2;
        field[nn] = '\0';
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
    check("fields", ifield > 40);
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
    same_fields(h1, h2);
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
