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
 * running one TOML input, with the same sizes and no failure, and a run
 * of a broken input failing as an input failure at the input stage.
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
