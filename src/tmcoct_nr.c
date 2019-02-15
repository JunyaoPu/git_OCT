#include <stdlib.h>
#include <stdio.h>
#include "tmcoct.h"

/***********************************************************
 *	Report error message to stderr, then exit the program
 *	with signal 1.
 ****/
void nrerror(char error_text[]) {
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}

/***********************************************************
 *	Allocate a matrix with row index from nrl to nrh 
 *	inclusive, and column index from ncl to nch
 *	inclusive.
 ****/
double **AllocMatrix(short nrl, short nrh,
                     short ncl, short nch) {
    short i, j;
    double **m;

    m = (double **) malloc((unsigned) (nrh - nrl + 1)
                           * sizeof(double *));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++) {
        m[i] = (double *) malloc((unsigned) (nch - ncl + 1)
                                 * sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in matrix()");
        m[i] -= ncl;
    }

    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++) m[i][j] = 0.0;
    return m;
}


double ***AllocHyperMatrix(short nrl, short nrh,
                           short ncl, short nch,
                           short nzl, short nzh) {
    short i, j, k;
    double ***m;
    int x = nrh - nrl;
    int y = nch - ncl;
    int z = nzh - nzl;

    m = (double ***) malloc(x * sizeof(double **));

    if (!m) {
        printf("\n\nError!  Not enough memory!\n\n");
        exit(1);
    }

//create second subscript
    for (i = 0; i < x; i++) {
        m[i] = (double **) malloc(y * sizeof(double *));

        if (!m[i]) {
            printf("\n\nError!  Not enough memory!\n\n");
            exit(1);
        }

        //create third subscript
        for (j = 0; j < y; j++) {
            m[i][j] = (double *) malloc(z * sizeof(double));

            if (!m[i][j]) {
                printf("\n\nError!  Not enough memory!\n\n");
                exit(1);
            }

        }

    }
    for (i = 0; i < x; i++)
        for (j = 0; j < y; j++)
            for (k = 0; k < z; k++)
                m[i][j][k] = 0;


    return m;
}

/***********************************************************
 *	Release the memory.
 ****/
void FreeVector(double *v, short nl, short nh) {
    free((char *) (v + nl));
}

void FreeMatrix(double **m, short nrl, short nrh,
                short ncl, short nch) {
    short i;

    for (i = nrh; i >= nrl; i--) free((char *) (m[i] + ncl));
    free((char *) (m + nrl));
}

/***********************************************************
 *	Release the memory.
 ****/
void FreeHyperMatrix(double ***m, short nrl, short nrh,
                     short ncl, short nch, short nzl, short nzh) {
    short i, j;

    for (i = nrl; i < nrh; i++) {

        for (j = ncl; j < nch; j++)
            free(m[i][j]);

        free(m[i]);
    }

    free(m);
}
