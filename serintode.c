
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "gmp.h"
#include "iml.h"

//Check things with Valgrind
//Compile: gcc -Wall serintode.c -o serintode -liml -lcblas -lgmp -lm

int main()
{
    long const MAX_DIGITS=1000L; /*This needs to be input*/
    long const NUM_COEFFS=31L; /*This needs to be input. If changing type, change n declaration and n equating twice below*/ 
    long const MAX_ODE_ORDER=2L; /*This needs to be input or automated*/
    long const MAX_POLY_ORDER=3L; /*This needs to be input or automated*/
    long const COLUMNS=(MAX_ODE_ORDER+1)*(MAX_POLY_ORDER+1);
    long const ROWS=COLUMNS+5L; //5 checks
    long i,j,k;
    long nulldim;
    char input_string[MAX_DIGITS+1L];
    mpz_t *S, *M, *N, temp, temp2,coeff;
    FILE* fid=NULL;
    char *fgcheck;
    
    mpz_inits(temp,temp2,coeff,NULL);

    S = (mpz_t*) malloc(NUM_COEFFS*sizeof(mpz_t));
    
    fid = fopen("tests/catalan.txt","r");
    if (fid==NULL)
    {
        printf("\nERROR: Could not open input file. %s\n",strerror(errno));
        return 1;
    }
    
    for (i=0L;i<NUM_COEFFS;i++)
    {
        fgcheck = fgets(input_string,MAX_DIGITS,fid);
        if (fgcheck==NULL)
        {
            if (feof(fid))
            {
                printf("End of file reached before total number %ld of coefficients read.\n",NUM_COEFFS);
            }
            else if (ferror(fid))
            {
                printf("Error while reading the file. %s\n",strerror(errno));
            }
            for (j=0L;j<i+1L;j++)
            {
                mpz_clear(S[j]);
            }
            fclose(fid);
            return 1;
        }    
        input_string[strcspn(input_string,"\r\n")]='\0'; /*Remove trailing newline characters, May not be necessary*/
        mpz_init_set_str(S[i],input_string,10);
        //mpz_out_str(NULL,10,M[i]);
        //printf("\n");
    }
    
    M = (mpz_t*) malloc(ROWS*COLUMNS*sizeof(mpz_t));
    
    for (i=0L;i<MAX_POLY_ORDER+1;i++)
    {
        for (j=0;j<MAX_ODE_ORDER+1;j++)
        {
            for (k=i;k<ROWS;k++)
            {
                mpz_fac_ui(temp,j+k-i);
                mpz_fac_ui(temp2,k-i);
                mpz_divexact(coeff,temp,temp2); //coeff = (j+k-i)!/(k-i)!
                mpz_mul(temp,coeff,S[j+k-i]);
                mpz_init_set(M[k*COLUMNS+i*MAX_POLY_ORDER+j],temp);
            }
        }
    }
    
    /*
    printf("Input Matrix M:\n");
    for (i = 0; i < ROWS; i++)
    {
        for (j = 0; j < COLUMNS; j++)
        {
            gmp_fprintf (stdout, "  %Zd", M[i * COLUMNS + j]);
        }
        fprintf (stdout, "\n");
    }*/
    
    
    nulldim = nullspaceMP (ROWS, COLUMNS, M, &N);
    fprintf (stdout, "Dimension of nullspace: ");
    fprintf (stdout, " %ld\n", nulldim);
    for (i = 0L; i < COLUMNS; i++)
    {
        for (j = 0L; j < nulldim; j++)
        {
            gmp_fprintf (stdout, "  %Zd", N[i * nulldim + j]);
        }
        fprintf (stdout, "\n"); 
    }
    
    
    /*Clearing mpz variables and closing file*/
    for (i=0L;i<NUM_COEFFS;i++)
    {
        mpz_clear(S[i]);
    }
    for (i=0L;i<ROWS*COLUMNS;i++)
    {
        mpz_clear(M[i]);
    }
    for (i=0L;i<COLUMNS*nulldim;i++)
    {
        mpz_clear(N[i]);
    }
    mpz_clears(temp,temp2,coeff,NULL);
    free(M);
    free(N);
    fclose(fid);
    return 0;
}
