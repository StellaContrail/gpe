// COMPILE: gcc introduction_clang.c -lsilo -lm

#include <silo.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    DBfile *dbfile = NULL;
    /* Open the Silo file */
    dbfile = DBCreate("introduction_clang.silo", DB_CLOBBER, DB_LOCAL, "Test", DB_PDB);
    if(dbfile == NULL)
    {
        fprintf(stderr, "Could not create Silo file!\n");
        return -1;
    }
    /* Close the Silo file. */
    DBClose(dbfile);
    return 0;
}