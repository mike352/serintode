
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "iml.h"


//Compile: gcc -Wall serintode_nonlin2.c -o serintode_nonlin2.o -liml -lcblas -lgmp -lm
//Output file sum: gives the sequence name of found solutions, the number of coefficients, the ODE order, the largest polynomial order, the number of free variables
//Output file eqs: gives the sequence name of found solutions, plus the Maple input for the ODE

int main()
{
    char *finname = "tests/6-colorings.txt"; /*File name of data*/
    long const NUM_CHECKS=10L; /*Should be greater than 0*/
    long const MIN_ODE_ORDER=1L; 
    long const MAX_COEFFS=2000; /*Should be checked for very large sequences*/
    long const MAX_LINE_LENGTH=100000L; 
    //char foutsumname[64]; /*Output summary file name*/
    char fouteqsname[64]; /*Output equations file name*/
    long NUM_COEFFS=0L;
    long MAX_ODE_ORDER=0L; 
    long ODE_ORDER=0L;
    long MAX_POLY_ORDER=0L;
    long MAX_FOUND_ORDER=0L;
    long COLUMNS=0L, ROWS=0L;
    long i,j,k,l,m,n,numterms,maxnumterms;
    long nulldim=0L;
    char input_string[MAX_LINE_LENGTH+1L];
    mpz_t *I, **S, *M, *N, temp, temp2,coeff;
    FILE *fin=NULL, *fouteqs=NULL; //, *foutsum=NULL
    char *fgcheck, nulldimflag=0;
    long **orderexp;
    
    mpz_inits(temp,temp2,coeff,NULL);

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
    
    I = (mpz_t *) malloc(MAX_COEFFS*sizeof(mpz_t));
    M = (mpz_t *) malloc(MAX_COEFFS*MAX_COEFFS*sizeof(mpz_t)); //Maximum number needed below
    
    if((I == NULL) || (M == NULL))
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
                free(I);
                fclose(fin);
                exit(EXIT_FAILURE);
            }
        }   
        else
        {
            mpz_init_set_str(I[NUM_COEFFS],input_string,10);
            NUM_COEFFS++;
        }
    }
    fclose(fin);
    
    MAX_ODE_ORDER=floor((-5L+sqrt(17L+8L*(NUM_COEFFS-NUM_CHECKS)))/2L);     
    maxnumterms=(MAX_ODE_ORDER+1L)*(MAX_ODE_ORDER+2L)/2L;
    
    orderexp = (long **) malloc(maxnumterms*sizeof(long *));
    for (i=0L;i<maxnumterms;i++)
    {
        orderexp[i] = (long *) malloc(2L*sizeof(long));
    }
    
    //MAX_ODE_ORDER=floor((NUM_COEFFS-NUM_CHECKS-2L)/3L);
    for (n=MIN_ODE_ORDER;n<MAX_ODE_ORDER+1L;n++)
    {
        ODE_ORDER=n;
        numterms=(ODE_ORDER+1L)*(ODE_ORDER+2L)/2L;
        S = (mpz_t **) malloc(numterms*sizeof(mpz_t *));
        for (i=0L;i<numterms;i++)
        {
            S[i] = (mpz_t *) malloc(NUM_COEFFS*sizeof(mpz_t));
        }
        
        for (i=0L;i<numterms;i++)
        {
            for (j=0L;j<NUM_COEFFS;j++)
            {//printf("%ld %ld\n",i,j);
                mpz_init(S[i][j]);
            }
        }
        
        for (i=0L;i<NUM_COEFFS;i++)
        {
            mpz_set(S[0L][i],I[i]);
        }
        orderexp[0][0]=0;
        orderexp[0][1]=0;
        for (i=1L;i<numterms;i++)
        {
            for (j=0L;j<NUM_COEFFS-i;j++)
            {
                mpz_mul_ui(temp,S[i-1L][j+1L],j+1L);
                mpz_set(S[i][j],temp);
            }
            orderexp[i][0]=i;
            orderexp[i][1]=0;
        }
        m=0L;
        for (i=0L;i<ODE_ORDER;i++)
        {
            for (j=i;j<ODE_ORDER;j++)
            {
                m++;
                for (k=0L;k<NUM_COEFFS-1L;k++)
                {
                    mpz_set_ui(temp,0L);
                    for (l=0L;l<k+1L;l++)
                    {
                        mpz_mul(temp2,S[i][l],S[j][k-l]);
                        mpz_add(temp,temp,temp2);
                    }
                    mpz_set(S[ODE_ORDER+m][k],temp);
                } 
                orderexp[ODE_ORDER+m][0]=i;
                orderexp[ODE_ORDER+m][1]=j;
            }
        }
        
        MAX_POLY_ORDER=floor((NUM_COEFFS-NUM_CHECKS-ODE_ORDER)/numterms)-1L;
        COLUMNS=numterms*(MAX_POLY_ORDER+1L);
        ROWS=COLUMNS+NUM_CHECKS; 
        printf("Checking %s for nonlinear ODE of ""order"" %ld with polynom coeffs of order %ld and %ld checks\n",finname,ODE_ORDER,MAX_POLY_ORDER,NUM_CHECKS);
        
        for (i=0L;i<COLUMNS;i++)
        {
            for (j=0L;j<ROWS;j++)
            {
                mpz_init(M[j*COLUMNS+i]);
            }
        }
        
        for (i=0L;i<MAX_POLY_ORDER+1L;i++)
        {
            for (j=i;j<ROWS;j++)
            {
                for (k=0L;k<numterms;k++)
                {
                    mpz_set(M[j*COLUMNS+i*numterms+k],S[k][j-i]);
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
            //Check that the order is not 0
            mpz_set_ui(temp,0L);
            for (i=1L;i<numterms;i++)
            {
                for (j=0L;j<MAX_POLY_ORDER+1L;j++)
                {
                    mpz_abs(temp2,N[(i+j*numterms)*nulldim]);
                    mpz_add(temp,temp,temp2);
                }
            }
            if (mpz_get_ui(temp)==0L)
            {
                //printf("Spurious order 0 equation. ");
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
    
        for (i=0L;i<numterms;i++)
        {
            for (j=0L;j<NUM_COEFFS;j++)
            {
                mpz_clear(S[i][j]);
            }
        }
        for (i=0L;i<numterms;i++)
        {
            free(S[i]);
        }
        free(S);
    }
    
    
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
        
        sprintf(fouteqsname,"%s_nonlinsol_%ld-checks.txt",finname,NUM_CHECKS);
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
            exit(EXIT_FAILURE);
        }
        setvbuf(fouteqs,NULL,_IOLBF,32);
        //fprintf(foutsum,"%s: %ld, %ld, ",finname,NUM_COEFFS,ODE_ORDER);
        fprintf(fouteqs,"ODE%s := ",finname);
        for (i=0L;i<numterms;i++)
        {
            mpz_set_ui(temp,0L);
            for (j=0L;j<MAX_POLY_ORDER+1L;j++) 
            {
                mpz_abs(temp2,N[(i+j*numterms)*nulldim]);
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
                    if (mpz_cmp_ui(N[(i+k*numterms)*nulldim],0L)!=0L)
                    {
                        if ((k>0L)&&(mpz_cmp_ui(N[(i+k*numterms)*nulldim],0L)>0L))
                        {
                            fprintf(fouteqs,"+");
                            printf("+");
                        }
                        gmp_fprintf (fouteqs, "%Zd", N[(i+k*numterms)*nulldim]);
                        gmp_fprintf (stdout, "%Zd", N[(i+k*numterms)*nulldim]);
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
                if (i==0L)
                {
                    fprintf(fouteqs,"*y(x)");
                    printf("*y");
                }
                else if (i==1L)
                {
                    fprintf(fouteqs,"*diff(y(x),x)");
                    printf("*Dx");
                }
                else if (i<ODE_ORDER+1L)
                {
                    fprintf(fouteqs,"*diff(y(x),x$%ld)",i);
                    printf("*Dx^%ld",i);
                }
                else if (orderexp[i][0]==0)
                {
                    if (orderexp[i][1]==0)
                    {
                        fprintf(fouteqs,"*(y(x))^2");
                        printf("*y^2");
                    }
                    else
                    {
                        fprintf(fouteqs,"*(y(x))*(diff(y(x),x$%ld))",orderexp[i][1]);
                        printf("*y*Dx^%ld",orderexp[i][1]);
                    }
                }
                else
                {
                    if (orderexp[i][1]==orderexp[i][0])
                    {
                        fprintf(fouteqs,"*(diff(y(x),x$%ld))^2",orderexp[i][1]);
                        printf("*(Dx^%ld)^2",orderexp[i][1]);
                    }
                    else
                    {
                        fprintf(fouteqs,"*(diff(y(x),x$%ld))*(diff(y(x),x$%ld))",orderexp[i][0],orderexp[i][1]);
                        printf("*Dx^%ld*Dx%ld",orderexp[i][0],orderexp[i][1]);
                    }
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
                
    for (i=0L;i<maxnumterms;i++)
    {
        free(orderexp[i]);
    }
    free(orderexp);
    mpz_clears(temp,temp2,coeff,NULL);
    //fclose(foutsum);
    return 0;
}
