
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gmp.h"
#include "iml.h"


//Compile: gcc -Wall serintode_iml_auto.c -o serintode_iml_auto.o -liml -lcblas -lgmp -lm
//Output file sum: gives the sequence name of found solutions, the number of coefficients, the ODE order, the largest polynomial order, the number of free variables
//Output file eqs: gives the sequence name of found solutions, plus the Maple input for the ODE

int main()
{
    time_t start,end;
    char *finname = "tests/new_Chi3_w.ser"; /*File name of data*/
    long const NUM_CHECKS=10L; /*Should be greater than 0*/
    long const MIN_ODE_ORDER=1L; 
    long const MAX_COEFFS=400; /*Should be checked for very large sequences*/
    long const MAX_LINE_LENGTH=100000L; 
    //char foutsumname[64]; /*Output summary file name*/
    char fouteqsname[64]; /*Output equations file name*/
    long NUM_COEFFS=0L;
    long MAX_ODE_ORDER=0L; 
    long ODE_ORDER=0L;
    long MAX_POLY_ORDER=0L;
    long MAX_FOUND_ORDER=0L;
    long COLUMNS=0L, ROWS=0L;
    long i,j,k,n,nonzeroterms;
    long nulldim=0L;
    char input_string[MAX_LINE_LENGTH+1L];
    mpz_t *S, *M, *N, temp, temp2,coeff;
    FILE *fin=NULL, *fouteqs=NULL; //, *foutsum=NULL
    char *fgcheck, nulldimflag=0;
    
    time(&start);
    
    mpz_inits(temp,temp2,coeff,NULL);
    
    /*
    for (i=0L;i<MAX_COEFFS;i++)
    {
        mpz_init(S[i]);
    }
    for (i=0L;i<MAX_COEFFS*MAX_COEFFS;i++)
    {
        mpz_init(M[i]);
    }
    */
    
    fin = fopen(finname,"r");
    if (fin==NULL)
    {
        printf("\nError: Could not open input file %s. %s\n",finname,strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    /*
    sprintf(foutsumname,"%s_summary_%ld-checks.txt",finname,NUM_CHECKS);
    foutsum = fopen(foutsumname,"w");
    if (foutsum==NULL)
    {
        printf("\nError: Could not open summary output file %s. %s\n",foutsumname,strerror(errno));
        fclose(fin);
        exit(EXIT_FAILURE);
    }
    setvbuf(foutsum,NULL,_IOLBF,32);
    */
    
    S = (mpz_t*) malloc(MAX_COEFFS*sizeof(mpz_t));
    M = (mpz_t*) malloc(MAX_COEFFS*MAX_COEFFS*sizeof(mpz_t)); //Maximum number needed below
    
    if((S == NULL) || (M == NULL))
    {
        fprintf(stderr, "No memory left for allocating with malloc for input matrices. %s",strerror(errno));
        fclose(fin);
        //fclose(foutsum);
        //fclose(fouteqs);
        exit(EXIT_FAILURE);
    }
    
    NUM_COEFFS=0L;
    while (NUM_COEFFS<MAX_COEFFS)
    {
        fgcheck = fgets(input_string,MAX_LINE_LENGTH,fin);
        if (fgcheck==NULL)
        {
            if (feof(fin))
            {
                break;
            }
            else if (ferror(fin))
            {
                printf("Error while reading the file. %s\n",strerror(errno));
                free(S);
                fclose(fin);
                exit(EXIT_FAILURE);
            }
        }   
        else
        {
            mpz_init_set_str(S[NUM_COEFFS],input_string,10);
            NUM_COEFFS++;
        }
    }
    fclose(fin);
        
    /*    
    printf("Input S Matrix:\n");
    for (i=0L;i<NUM_COEFFS;i++)
    {
        gmp_fprintf (stdout, "%Zd ", S[i]);
    }
    printf("\n\n");
    */
    
    printf("File %s:\n",finname);
    MAX_ODE_ORDER=floor((NUM_COEFFS-NUM_CHECKS)/2L)-1L;
    for (n=MIN_ODE_ORDER;n<MAX_ODE_ORDER+1L;n++)
    {
        ODE_ORDER=n; //printf("Starting ODE order %ld\n",ODE_ORDER);
        MAX_POLY_ORDER=floor((NUM_COEFFS-NUM_CHECKS)/(ODE_ORDER+1L))-2L;
        if (MAX_POLY_ORDER == 0L)
        {
            break;
        }
        COLUMNS=(ODE_ORDER+1L)*(MAX_POLY_ORDER+1L);
        ROWS=COLUMNS+NUM_CHECKS; 
        printf("Checking for ODE order %ld with polynom coeffs of order %ld and %ld checks\n",ODE_ORDER,MAX_POLY_ORDER,NUM_CHECKS);
        
        for (i=0L;i<COLUMNS;i++)
        {
            for (j=0L;j<ROWS;j++)
            {
                mpz_init(M[j*COLUMNS+i]);
            }
        }
        
        //Null vector has the form: 0th order poly coeffs of first term up to ODE order term, order 1 poly coeffs of all terms, order 2 poly coeffs of all terms, etc
        //Matrix rows correspond to each successive order that needs to be annihilated
        //Matrix columns are the result of the derivatives and poly coeffs acting on input series to give a particular final order
        for (i=0L;i<MAX_POLY_ORDER+1L;i++)
        {
            for (j=0;j<ODE_ORDER+1L;j++)
            {
                for (k=i;k<ROWS;k++)
                {
                    mpz_fac_ui(temp,j+k-i);
                    mpz_fac_ui(temp2,k-i);
                    mpz_divexact(coeff,temp,temp2); //coeff = (j+k-i)!/(k-i)!
                    mpz_mul(temp,coeff,S[j+k-i]); //printf("Check i=%ld, j=%ld, k=%ld, index=%ld, index=%ld, max=%ld  ",i,j,k,j+k-i,k*COLUMNS+i*(ODE_ORDER+1)+j,ROWS*COLUMNS);
                    //mpz_out_str(NULL,10,coeff); printf("\n");
                    mpz_set(M[k*COLUMNS+i*(ODE_ORDER+1L)+j],temp);
                }
            }
        }
        
        /*
        printf("Input %ld x %ld M matrix\n",ROWS,COLUMNS);
        for (i=0L;i<MAX_POLY_ORDER+1L;i++)
        {
            for (j=0;j<ODE_ORDER+1L;j++)
            {
                for (k=i;k<ROWS;k++)
                {
                    gmp_fprintf (stdout, "  %Zd", M[k*COLUMNS+i*(ODE_ORDER+1)+j]);
                }
            }
        }
        printf("\n\n");
        */
        
        /*
        printf("Input %ld x %ld M matrix\n",ROWS,COLUMNS);
        for (i=0L;i<ROWS;i++)
        {
            for (j=0L;j<COLUMNS;j++)
            {
                gmp_fprintf (stdout, "  %Zd", M[i*COLUMNS+j]);
            }
            printf("\n");
        }
        */
        
        //nulldim = nullspaceMP(ROWS,COLUMNS,M,&N);
        nulldim = kernelMP(ROWS,COLUMNS,M,&N,1L);
        
        for (i=0L;i<COLUMNS;i++)
        {
            for (j=0;j<ROWS;j++)
            {
                mpz_clear(M[j*COLUMNS+i]);
            }
        }
        
        /*
        printf("Output N matrix\n");
        for (i=0L;i<COLUMNS*nulldim;i++)
        {
            gmp_fprintf (stdout, "  %Zd", N[i]);
        }
        printf("\n\n");
        */
        
        mpz_set_ui(temp,0L);
        if (nulldim>0L)
        {
            //Checking whether kernel solution is spurious
            for (i = 0L; i < COLUMNS; i++)
            {
                for (j = 0L; j < nulldim; j++)
                {
                    mpz_abs(temp2,N[i * nulldim + j]);
                    mpz_add(temp,temp,temp2);
                }
            }
            if (mpz_get_ui(temp)!=nulldim)
            {
                nulldimflag=1;
            }
        }
        
        if (nulldimflag==1)
        {
            nonzeroterms=0L;
            for (i=0L;i<(ODE_ORDER+1L);i++)
            {
                mpz_set_ui(temp,0L);
                for (j=0L;j<MAX_POLY_ORDER+1L;j++)
                {
                    mpz_abs(temp2,N[(i+j*(ODE_ORDER+1L))*nulldim]);
                    mpz_add(temp,temp,temp2);
                }
                if (mpz_get_ui(temp)!=0L)
                {
                    nonzeroterms++;
                }
            }
            if (nonzeroterms<2L)
            {
                printf("Spurious equation with none or only 1 non-zero term.\n");
                nulldimflag=0;
                for (i=0L;i<COLUMNS*nulldim;i++)
                {
                    mpz_clear(N[i]);
                }
                free(N);
            }
            else
            {
                break;
            }
        }
        else
        {
            for (i=0L;i<COLUMNS*nulldim;i++)
            {
                mpz_clear(N[i]);
            }
            free(N);
        }
        
    }
    
    
    for (i=0L;i<NUM_COEFFS;i++)
    {
        mpz_clear(S[i]);
    }
    free(S);
    free(M);
    
    
    if (nulldimflag==1)
    {
        /*
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
        printf("\n\n");
        for (i = 0L; i<ODE_ORDER+1L; i++)
        {
            for (j = 0L; j < MAX_POLY_ORDER+1L; j++)
            {
                gmp_fprintf (stdout, "  %Zd", N[(i+j*(ODE_ORDER+1L))*nulldim]);
            }
            fprintf (stdout, "\n"); 
        }
        */
        printf("****Found a solution!****\n");
        
        sprintf(fouteqsname,"%s_solution_%ld-checks.txt",finname,NUM_CHECKS);
        fouteqs = fopen(fouteqsname,"w");
        if (fouteqs==NULL)
        {
            printf("\nError: Could not open equations output file %s. %s\n",fouteqsname,strerror(errno));
            fclose(fin);
            //fclose(foutsum);
            for (i=0L;i<COLUMNS*nulldim;i++)
            {
                mpz_clear(N[i]);
            }
            free(N);
            free(M);
            exit(EXIT_FAILURE);
        }
        setvbuf(fouteqs,NULL,_IOLBF,32);
        //fprintf(foutsum,"%s: %ld, %ld, ",finname,NUM_COEFFS,ODE_ORDER);
        fprintf(fouteqs,"ODE%s := ",finname);
        for (i=0L;i<ODE_ORDER+1L;i++)
        {
            mpz_set_ui(temp,0L);
            for (j=0L;j<MAX_POLY_ORDER+1L;j++) //Start from 1 to remove order 0 ODEs
            {
                mpz_abs(temp2,N[(i+j*(ODE_ORDER+1L))*nulldim]);
                mpz_add(temp,temp,temp2);
            }
            if (mpz_cmp_ui(temp,0L)>0L)
            {
                if (i>0L)
                {
                    fprintf(fouteqs,"+");
                    printf("+");
                }
                fprintf(fouteqs,"(");
                printf("(");
                for (k=0L;k<MAX_POLY_ORDER+1L;k++)
                {
                    if (mpz_cmp_ui(N[(i+k*(ODE_ORDER+1L))*nulldim],0L)!=0L)
                    {
                        if ((k>0L)&&(mpz_cmp_ui(N[(i+k*(ODE_ORDER+1L))*nulldim],0L)>0L))
                        {
                            fprintf(fouteqs,"+");
                            printf("+");
                        }
                        gmp_fprintf (fouteqs, "%Zd", N[(i+k*(ODE_ORDER+1L))*nulldim]);
                        gmp_fprintf (stdout, "%Zd", N[(i+k*(ODE_ORDER+1L))*nulldim]);
                        if (k>0L)
                        {
                            fprintf(fouteqs,"*x^%ld",k);
                            printf("*x^%ld",k);
                        }
                        if (MAX_FOUND_ORDER<k)
                        {
                            MAX_FOUND_ORDER=k;
                        }
                    }
                }
                fprintf(fouteqs,")");
                printf(")");
                if (i>0L)
                {
                    fprintf(fouteqs,"*diff(y(x),x$%ld)",i);
                    printf("*Dx^%ld",i);
                }
                else
                {
                    fprintf(fouteqs,"*y(x)");
                    printf("*y(x)");
                }
            }
        }
        fprintf(fouteqs,":\n");
        printf("=0\n");
        
        //fprintf(foutsum,"%ld, %ld\n",MAX_FOUND_ORDER,NUM_COEFFS-NUM_CHECKS-(ODE_ORDER+1L)*(MAX_FOUND_ORDER+1L));
        
        for (i=0L;i<COLUMNS*nulldim;i++)
        {
            mpz_clear(N[i]);
        }
        free(N);
        fclose(fouteqs);
    }
    else
    {
        if (ODE_ORDER>0L)
        {
            printf("Couldn't find any solution\n");
        }
        else
        {
            printf("Series too short %ld checks\n",NUM_CHECKS);
        }
    }
    
    time(&end);
    printf("\nEllapsed time %.fs\n",difftime(end,start));

    mpz_clears(temp,temp2,coeff,NULL);
    exit(EXIT_SUCCESS);
}
