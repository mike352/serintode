
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "gmp.h"
#include "iml.h"

//Check things with Valgrind
//Compile: gcc -Wall imltest.c -o imltest -liml -lcblas -lgmp -lm

int main()
{
    long const MAX_DIGITS=300L; /*This needs to be input*/
    long const NUM_INPUT=100L; /*This needs to be input. If changing type, change n declaration and n equating twice below*/ 
    long const ROWS=2L;
    long const COLUMNS=2L;
    long n,m,i,j;
    long s;
    char input_string[MAX_DIGITS+1L];
    mpz_t *v, *r;
    FILE* fid=NULL;
    char *fgcheck;
    

    v = (mpz_t*) malloc(NUM_INPUT*sizeof(mpz_t));
    
    fid = fopen("kolakoski.txt","r");
    if (fid==NULL)
    {
        printf("\nERROR: Could not open input file. %s\n",strerror(errno));
        return 1;
    }
    
    for (n=0L;n<NUM_INPUT;n++)
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
            for (m=0L;m<n+1L;m++)
            {
                mpz_clear(v[m]);
            }
            fclose(fid);
            return 1;
        }    
        input_string[strcspn(input_string,"\r\n")]='\0'; /*Remove trailing newline characters, May not be necessary*/
        mpz_init_set_str(v[n],input_string,10);
        //mpz_out_str(NULL,10,v[n]);
        //printf("\n");
    }
    
    n=ROWS; m=COLUMNS;
    mpz_set_si(v[0],2); //1
    mpz_set_si(v[1],1); //2
    mpz_set_si(v[2],-4); //3
    mpz_set_si(v[3],-2); //6
    s = nullspaceMP (n, m, v, &r);
    fprintf (stdout, "Dimension of nullspace: ");
    fprintf (stdout, " %ld\n", s);
    for (i = 0; i < m; i++)
        {
        for (j = 0; j < s; j++)
            gmp_fprintf (stdout, "  %Zd", r[i * s + j]);
        fprintf (stdout, "\n");
        }
    
    //Factorial test
    mpz_fac_ui(v[2],100);
    //mpz_out_str(NULL,10,v[2]);
    //printf("\n");
    
    /*Clearing mpz variables and closing file*/
    for (n=0L;n<NUM_INPUT;n++)
    {
        mpz_clear(v[n]);
    }
    for (n=0L;n<ROWS*s;n++)
    {
        mpz_clear(r[n]);
    }
    free(v);
    free(r);
    fclose(fid);
    return 0;
}
