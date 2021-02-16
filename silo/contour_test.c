// gcc contour_test.c -lsilo -lm

#include <silo.h>
#include <stdio.h>
#include <stdlib.h>

void write_grids(DBfile *dbfile)
{
    const int N = 5;
    const int dims[] = {N, N, N};
    const int ndims = 3;
    float x[N];
    float y[N];
    float z[N];
    const float dh = 1;
    for (int i = 0; i < N; i++) {
        x[i] = dh * i;
        y[i] = dh * i;
        z[i] = dh * i;
    }
    float *coords[] = {x, y, z};
    DBPutQuadmesh(dbfile, "grid_mesh", NULL, coords, dims, ndims, DB_FLOAT, DB_COLLINEAR, NULL);
}

void write_contour(DBfile *dbfile)
{
    const int N = 5;
    const int dims[] = {N, N, N};
    const int ndims = 3;
    float nodal[N * N * N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                nodal[i + j * N + k * N * N] = i + j + k;
            }
        }
    }
    DBPutQuadvar1(dbfile, "nodal", "grid_mesh", nodal, dims, ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
}

int main(int argc, char *argv[])
{
    DBfile *dbfile = NULL;
    /* Open the Silo file */
    dbfile = DBCreate("contour.silo", DB_CLOBBER, DB_LOCAL, "Contour Plot Test", DB_PDB);
    if(dbfile == NULL)
    {
        fprintf(stderr, "Could not create Silo file!\n");
        return -1;
    }

    /* Add other Silo calls here. */
    write_grids(dbfile);
    write_contour(dbfile);

    /* Close the Silo file. */
    DBClose(dbfile);
    return 0;
}