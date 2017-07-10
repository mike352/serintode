
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
    long const MAX_DIGITS=300L; /*This needs to be input*/
    long const NUM_INPUT=100L; /*This needs to be input. If changing type, change n declaration and n equating twice below*/ 
    long const MAX_ODE_ORDER=2L; /*This needs to be input or automated*/
    long const MAX_POLY_ORDER=4L; /*This needs to be input or automated*/
    long const ROWS=2L;
    long const COLUMNS=2L;
    long i,j,r=ROWS,c=COLUMNS;
    long nulldim;
    char input_string[MAX_DIGITS+1L];
    mpz_t *M, *N;
    FILE* fid=NULL;
    char *fgcheck;
    

    M = (mpz_t*) malloc(NUM_INPUT*sizeof(mpz_t));
    
    fid = fopen("kolakoski.txt","r");
    if (fid==NULL)
    {
        printf("\nERROR: Could not open input file. %s\n",strerror(errno));
        return 1;
    }
    
    for (i=0L;i<NUM_INPUT;i++)
    {
        fgcheck = fgets(input_string,MAX_DIGITS,fid);
        if (fgcheck==NULL)
        {
            if (feof(fid))
            {
                printf("End of file reached before total number %ld of coefficients read.\n",NUM_INPUT);
            }
            else if (ferror(fid))
            {
                printf("Error while reading the file. %s\n",strerror(errno));
            }
            for (j=0L;j<i+1L;j++)
            {
                mpz_clear(M[j]);
            }
            fclose(fid);
            return 1;
        }    
        input_string[strcspn(input_string,"\r\n")]='\0'; /*Remove trailing newline characters, May not be necessary*/
        mpz_init_set_str(M[i],input_string,10);
        //mpz_out_str(NULL,10,M[i]);
        //printf("\n");
    }
    
    mpz_set_si(M[0],2); //1
    mpz_set_si(M[1],1); //2
    mpz_set_si(M[2],-4); //3
    mpz_set_si(M[3],-2); //6
    nulldim = nullspaceMP (r, c, M, &N);
    fprintf (stdout, "Dimension of nullspace: ");
    fprintf (stdout, " %ld\n", nulldim);
    for (i = 0; i < c; i++)
        {
        for (j = 0; j < nulldim; j++)
            gmp_fprintf (stdout, "  %Zd", N[i * nulldim + j]);
        fprintf (stdout, "\n");
        }
    
    //Factorial test
    mpz_fac_ui(M[2],100);
    //mpz_out_str(NULL,10,M[2]);
    //printf("\n");
    
    /*Clearing mpz variables and closing file*/
    for (i=0L;i<NUM_INPUT;i++)
    {
        mpz_clear(M[i]);
    }
    for (i=0L;i<r*nulldim;i++)
    {
        mpz_clear(N[i]);
    }
    free(M);
    free(N);
    fclose(fid);
    return 0;
}
