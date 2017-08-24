
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "gmp.h"
#include "iml.h"


//Compile: gcc -Wall serintode_iml_nonlin.c -o nonlin.o -liml -lcblas -lgmp -lm
//Output file sum: gives the sequence name of found solutions, the number of coefficients, the ODE order, the largest polynomial order, the number of free variables
//Output file eqs: gives the sequence name of found solutions, plus the Maple input for the ODE


//Thanks to Richard Brent from Newcastle University for this function
void combs(long **r, long *s, long k, long p, long q, long *arrindex)
{
    long j, m;
    if (k <= 0)
    {
        for (m=0; m < q; m++)
    {
        r[*arrindex][m]=s[m];
    }
        (*arrindex)++;  
    }
    else
    {
        for (j=p; j>=0; j--)
        {
            s[q] = j;
            combs(r,s,k-1,p-j,q+1,arrindex);
        }
    }
    return;  
}

int main()
{
    char *finname = "tests/3-colorings.txt"; /*File name of data*/
    long const NUM_CHECKS=10L; /*Should be greater than 0*/
    long const MIN_ODE_ORDER=1L; 
    long const MIN_DEPTH=1L;
    long const MAX_DEPTH=10L; //Choosing 1 equals linear
    long const MAX_COEFFS=40; /*Should be checked for very large sequences*/
    long const MAX_LINE_LENGTH=100000L; 
    //char foutsumname[64]; /*Output summary file name*/
    char fouteqsname[64]; /*Output equations file name*/
    long NUM_COEFFS=0L;
    long MAX_ODE_ORDER=100L; 
    long ODE_ORDER=0L;
    long MAX_POLY_ORDER=0L;
    long MAX_FOUND_ORDER=0L;
    long COLUMNS=0L, ROWS=0L;
    long i,j,k,l,m,n,p,numterms=0L,maxnumterms=0L,ordermaxnumterms=0L,first,MAX_DEPTH_POSS,nonzeroterms,ordersused,termsused, MAX_FOUND_POLY_ORDER,firstterm, firstorder,mintermsused,MIN_MAX_FOUND_POLY_ORDER,MIN_MAX_FOUND_DEPTH, bestnulldim,finalorders, MAX_FOUND_DEPTH,MAX_FOUND_ODE_ORDER;
    long nulldim=0L;
    long **orderexp, *s, arrindex=0L;
    char input_string[MAX_LINE_LENGTH+1L];
    mpz_t *I, **D, **S, *M, *N, temp, temp2,coeff, *temparray;
    FILE *fin=NULL, *fouteqs=NULL; //, *foutsum=NULL
    char *fgcheck, nulldimflag=0;
    time_t start,end;
    
    time(&start);
    setvbuf(stdout,NULL,_IONBF,0);
    
    mpz_inits(temp,temp2,coeff,NULL);

    fin = fopen(finname,"r");
    if (fin==NULL)
    {
        printf("\nError: Could not open input file %s. %s\n",finname,strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    
    I = (mpz_t *) malloc(MAX_COEFFS*sizeof(mpz_t));
    M = (mpz_t *) malloc(MAX_COEFFS*MAX_COEFFS*sizeof(mpz_t)); //Maximum number needed below
    temparray = (mpz_t *) malloc(MAX_COEFFS*sizeof(mpz_t));
    
    if((I == NULL) || (M == NULL) || (temparray==NULL))
    {
        fprintf(stderr, "No memory left for allocating I, M, or temparray. %s",strerror(errno));
        fclose(fin);
        exit(EXIT_FAILURE);
    }
    
    for (i=0L;i<MAX_COEFFS;i++)
    {
        mpz_init(temparray[i]);
    }
    
    NUM_COEFFS=0L;
    nonzeroterms=0L;
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
                free(M);
                free(temparray);
                fclose(fin);
                exit(EXIT_FAILURE);
            }
        }   
        else
        {
            if (nonzeroterms==0L) //Remove initial zeros
            {
                mpz_set_str(temp,input_string,10);
                if (mpz_cmp_ui(temp,0)!=0)
                {
                    nonzeroterms++;
                    mpz_init_set(I[NUM_COEFFS],temp);
                    NUM_COEFFS++;
                }
            }
            else
            {
                mpz_init_set_str(I[NUM_COEFFS],input_string,10);
                NUM_COEFFS++;
            }
        }
    }
    fclose(fin); 
    printf("\nChecking %s:\n",finname);
    
    MAX_ODE_ORDER=NUM_COEFFS-NUM_CHECKS-1L;     
    maxnumterms=NUM_COEFFS-NUM_CHECKS;
    
    D = (mpz_t **) malloc((MAX_ODE_ORDER+1L)*sizeof(mpz_t *));
    if(D==NULL)
    {
        fprintf(stderr, "No memory left for allocating D matrix. %s",strerror(errno));
        free(I);
        free(M);
        free(temparray);
        exit(EXIT_FAILURE);
    }
    for (i=0L;i<(MAX_ODE_ORDER+1L);i++)
    {
        D[i] = (mpz_t *) calloc(NUM_COEFFS,sizeof(mpz_t));
        if(D[i]==NULL)
        {
            fprintf(stderr, "No memory left for allocating D array. %s",strerror(errno));
            for (j=0L;j<i;j++)
            {
                free(D[j]);
            }
            free(D);
            free(I);
            free(M);
            free(temparray);
            exit(EXIT_FAILURE);
        }
    }
    for (i=0L;i<(MAX_ODE_ORDER+1L);i++)
    {
        for (j=0L;j<NUM_COEFFS;j++)
        {
            mpz_init(D[i][j]);
        }
    }
    
    for (i=0L;i<NUM_COEFFS;i++)
    {
        mpz_set(D[0L][i],I[i]);
    }
    for (i=1L;i<MAX_ODE_ORDER+1L;i++)
    {
        for (j=0L;j<NUM_COEFFS-i;j++)
        {
            mpz_mul_ui(temp,D[i-1L][j+1L],j+1L);
            mpz_set(D[i][j],temp);
        }
    }
    
    /*
    for (i=0L;i<3L;i++)
    {
        printf("D[%ld]: ",i);
        for (j=0L;j<NUM_COEFFS-i;j++)
        {
            gmp_printf("%Zd ",D[i][j]);
        }
        printf("\n");
    }*/
    
    //MAX_ODE_ORDER=floor((double) (NUM_COEFFS-NUM_CHECKS-2L)/3L);
    for (n=MIN_ODE_ORDER;n<MAX_ODE_ORDER+1L;n++)
    {
        ODE_ORDER=n;
        for (MAX_DEPTH_POSS=1;MAX_DEPTH_POSS<=maxnumterms;MAX_DEPTH_POSS++)
        {
            //Computing the number (n+d)!/(n!d!). Need one less than this. n is order d is depth
            mpz_fac_ui(temp,ODE_ORDER+1L+MAX_DEPTH_POSS);
            mpz_fac_ui(temp2,ODE_ORDER+1L);
            mpz_divexact(temp,temp,temp2);
            mpz_fac_ui(temp2,MAX_DEPTH_POSS);
            mpz_divexact(temp,temp,temp2);
            if ((long) mpz_get_ui(temp)>maxnumterms)
            {
                MAX_DEPTH_POSS--;
                break;
            }
            else
            {
                ordermaxnumterms=mpz_get_ui(temp);
            }
        }
        if (MAX_DEPTH_POSS>MAX_DEPTH)
        {
            MAX_DEPTH_POSS=MAX_DEPTH;
        }
        
        S = (mpz_t **) malloc(ordermaxnumterms*sizeof(mpz_t *));
        if(S==NULL)
        {
            fprintf(stderr, "No memory left for allocating S matrix. %s",strerror(errno));
            for (j=0L;j<NUM_COEFFS;j++)
            {
                free(D[j]);
            }
            free(D);
            free(I);
            free(M);
            free(temparray);
            exit(EXIT_FAILURE);
        }
        for (i=0L;i<ordermaxnumterms;i++)
        {
            S[i] = (mpz_t *) calloc(NUM_COEFFS,sizeof(mpz_t));
            if(S[i]==NULL)
            {
                fprintf(stderr, "No memory left for allocating S array. %s",strerror(errno));
                for (j=0L;j<NUM_COEFFS;j++)
                {
                    free(D[j]);
                }
                free(D);
                for (j=0L;j<i;j++)
                {
                    free(S[j]);
                }
                free(S);
                free(I);
                free(M);
                free(temparray);
                exit(EXIT_FAILURE);
            }
        }
        
        for (i=0L;i<ordermaxnumterms;i++)
        {
            for (j=0L;j<NUM_COEFFS;j++)
            {
                mpz_init(S[i][j]);
            }
        }

        //depthmax=1; //Linear ODE
        for (p=MIN_DEPTH;p<MAX_DEPTH_POSS+1L;p++)
        {
            mpz_fac_ui(temp,ODE_ORDER+1L+p);
            mpz_fac_ui(temp2,ODE_ORDER+1L);
            mpz_divexact(temp,temp,temp2);
            mpz_fac_ui(temp2,p);
            mpz_divexact(temp,temp,temp2);
            numterms=mpz_get_ui(temp)-1L; //Equals (n+p)!/(n!p!) Need one less than this
            MAX_POLY_ORDER=floor((double) (NUM_COEFFS-NUM_CHECKS-ODE_ORDER)/numterms)-1L;
            if (MAX_POLY_ORDER<1L)
            {
                break;
            }
            
            s = (long *) calloc((ODE_ORDER+1),sizeof(long));
            orderexp = (long **) malloc((numterms+1L)*sizeof(long *));
            if ((orderexp==NULL) || (s==NULL))
            {
                fprintf(stderr, "No memory left for allocating orderexp matrix. %s",strerror(errno));
                for (j=0L;j<NUM_COEFFS;j++)
                {
                    free(D[j]);
                }
                free(D);
                for (j=0L;j<ordermaxnumterms;j++)
                {
                    free(S[j]);
                }
                free(S);
                free(I);
                free(M);
                free(temparray);
                exit(EXIT_FAILURE);
            }
            for (i=0;i<(numterms+1L);i++)
            {
                orderexp[i] = (long *) calloc((ODE_ORDER+1),sizeof(long));
                if(orderexp[i]==NULL)
                {
                    fprintf(stderr, "No memory left for allocating orderexp array. %s",strerror(errno));
                    for (j=0L;j<NUM_COEFFS;j++)
                    {
                        free(D[j]);
                    }
                    free(D);
                    for (j=0L;j<ordermaxnumterms;j++)
                    {
                        free(S[j]);
                    }
                    free(S);
                    for (j=0L;j<i;j++)
                    {
                        free(orderexp[j]);
                    }
                    free(s);
                    free(I);
                    free(M);
                    free(temparray);
                exit(EXIT_FAILURE);
                }
            }
            arrindex=0L;
            
            combs(orderexp,s,ODE_ORDER+1L,p,0,&arrindex); 
            
            /*
            for (i=0L;i<numterms;i++)
            {
                for (j=0L;j<ODE_ORDER+1L;j++)
                {
                    printf("i=%ld, j=%ld, orderxp[i][j]=%ld\n",i,j,orderexp[i][j]);
                }
            }*/

            for (i=0L;i<numterms;i++)
            {
                first=0L;
                for (j=0L;j<ODE_ORDER+1L;j++)
                {
                    if (orderexp[i][j]==0L)
                    {
                        continue;
                    }
                    else if (first==0L)
                    {
                        for (k=0L;k<NUM_COEFFS-j;k++)
                        {
                            mpz_set(S[i][k],D[j][k]);
                        }
                        first++;
                        for (m=1L;m<orderexp[i][j];m++)
                        {
                            for (k=0L;k<NUM_COEFFS-j;k++)
                            {
                                mpz_set_ui(temp,0L);
                                for (l=0L;l<k+1L;l++)
                                {
                                    mpz_mul(temp2,S[i][l],D[j][k-l]);
                                    mpz_add(temp,temp,temp2);
                                }
                                mpz_set(temparray[k],temp);
                            }
                            for (k=0L;k<NUM_COEFFS-j;k++)
                            {
                                mpz_set(S[i][k],temparray[k]);
                            }
                        }
                    }
                    else
                    {
                        for (m=0L;m<orderexp[i][j];m++)
                        {
                            for (k=0L;k<NUM_COEFFS;k++)
                            {
                                mpz_set_ui(temp,0L);
                                for (l=0L;l<k+1L;l++)
                                {
                                    mpz_mul(temp2,S[i][l],D[j][k-l]);
                                    mpz_add(temp,temp,temp2);
                                }
                                mpz_set(temparray[k],temp);
                            }
                            for (k=0L;k<NUM_COEFFS-j;k++)
                            {
                                mpz_set(S[i][k],temparray[k]);
                            }
                        }
                    }
                }
            }
            
            COLUMNS=numterms*(MAX_POLY_ORDER+1L);
            ROWS=COLUMNS+NUM_CHECKS; 
            printf("ODE of order %ld with %ld nonlin terms depth %ld, polynom coeffs of order %ld and %ld checks, size %ld x %ld\n",ODE_ORDER,numterms,p,MAX_POLY_ORDER,NUM_CHECKS, ROWS,COLUMNS);
            
            for (i=0L;i<COLUMNS;i++)
            {
                for (j=0L;j<ROWS;j++)
                {
                    mpz_init(M[j*COLUMNS+i]);
                }
            }
            
            //Null vector has the form: 0th order poly coeffs of first term up to numterms, order 1 poly coeffs of all terms, order 2 poly coeffs of all terms, etc
            //Matrix rows correspond to each successive order that needs to be annihilated
            //Matrix columns are the result of the derivatives and poly coeffs acting on input series to give a particular final order
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
            nulldimflag=0L;
            nulldim = kernelMP(ROWS,COLUMNS,M,&N,1L);
            printf("null dimension = %ld\n",nulldim);
            
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
                else
                {
                    //printf("Spurious null vectors.\n");
                }
            }
            
            if (nulldimflag==1)
                    {
                        mintermsused=numterms+1L;
                        MIN_MAX_FOUND_POLY_ORDER=MAX_POLY_ORDER+1L;
                        MIN_MAX_FOUND_DEPTH=MAX_DEPTH+1L;
                        bestnulldim=0L;
                        finalorders=numterms+1L;
                        MAX_FOUND_DEPTH=0L;
                        for (k=0L;k<nulldim;k++)
                        {
                            //Check that the number of orders is not 1
                            ordersused=0L;
                            termsused=0L;
                            nonzeroterms=0L;
                            MAX_FOUND_POLY_ORDER=0L;
                            MAX_FOUND_ODE_ORDER=0L;
                            for (i=0L;i<numterms;i++)
                            {
                                for (j=0L;j<MAX_POLY_ORDER+1L;j++)
                                {
                                    if (mpz_cmp_ui(N[(i+j*numterms)*nulldim+k],0L)!=0L)
                                    {
                                        nonzeroterms++;
                                        if (MAX_FOUND_POLY_ORDER<j)
                                        {
                                            MAX_FOUND_POLY_ORDER=j;
                                        }
                                    }
                                }
                                if (termsused<nonzeroterms)
                                {
                                    termsused=nonzeroterms;
                                    ordersused++; //Not really order, but numterms
                                    for (l=0L;l<ODE_ORDER+1L;l++)
                                    {
                                        if (orderexp[i][l]!=0)
                                        {
                                            if (MAX_FOUND_ODE_ORDER<l)
                                            {
                                                MAX_FOUND_ODE_ORDER=l;
                                            }
                                            if (MAX_FOUND_DEPTH<orderexp[i][l])
                                            {
                                                MAX_FOUND_DEPTH=orderexp[i][l];
                                            }
                                        }
                                    }
                                }
                            }
                            if (MAX_FOUND_ODE_ORDER==ODE_ORDER)
                            {
                                if ((finalorders>ordersused)&&(finalorders>1L))
                                {
                                    finalorders=ordersused;
                                    MIN_MAX_FOUND_DEPTH=MAX_FOUND_DEPTH;
                                    MIN_MAX_FOUND_POLY_ORDER=MAX_FOUND_POLY_ORDER;
                                    mintermsused=termsused;
                                    bestnulldim = k;
                                }
                                else if ((finalorders==ordersused)&&(finalorders>1L))
                                {
                                    if (MIN_MAX_FOUND_DEPTH>MAX_FOUND_DEPTH)
                                    {
                                        MIN_MAX_FOUND_DEPTH=MAX_FOUND_DEPTH;
                                        MIN_MAX_FOUND_POLY_ORDER=MAX_FOUND_POLY_ORDER;
                                        mintermsused=termsused;
                                        bestnulldim = k;
                                    }
                                    else if (MIN_MAX_FOUND_DEPTH==MAX_FOUND_DEPTH)
                                    {
                                        if (MIN_MAX_FOUND_POLY_ORDER>MAX_FOUND_POLY_ORDER)
                                        {
                                            MIN_MAX_FOUND_POLY_ORDER=MAX_FOUND_POLY_ORDER;
                                            mintermsused=termsused;
                                            bestnulldim = k;
                                        }
                                    }
                                }
                            }
                        } //End of k nulldim loop
                        if (finalorders<2L)
                        {
                            //printf("Spurious equation with %ld order term.\n",ordersused);
                            nulldimflag=0;
                        }
                        else if (MAX_FOUND_ODE_ORDER<ODE_ORDER)
                        {
                            //printf("Spurious solution came from taking a lot of derivatives, then found lower order ODE.\n");
                            nulldimflag=0;
                        }
                        else if (mintermsused<finalorders+1L)
                        {
                            if (MIN_MAX_FOUND_DEPTH<2L) //Won't allow constat coeffs for linear ODEs
                            {
                                //printf("\nWARNING: Polynomial coefficients only have one term each.\n\n");
                                nulldimflag=0;
                            }
                        }
                        
                        if (nulldimflag==1)
                        {
                            break;
                        }
                        else
                        {
                            for (i=0L;i<COLUMNS*nulldim;i++)
                            {
                                mpz_clear(N[i]);
                            }
                            free(N);
                            for (i=0L;i<(numterms+1L);i++)
                            {
                                free(orderexp[i]);
                            }
                            free(orderexp);
                            free(s);
                        }
                    }
            else
            {
                for (i=0L;i<COLUMNS*nulldim;i++)
                {
                    mpz_clear(N[i]);
                }
                free(N);
                
                for (i=0L;i<(numterms+1L);i++)
                {
                    free(orderexp[i]);
                }
                free(orderexp);
                free(s);
            }
            //Repeat with greater depth
        } 
        
        for (i=0L;i<ordermaxnumterms;i++)
        {
            for (j=0L;j<NUM_COEFFS;j++)
            {
                mpz_clear(S[i][j]);
            }
        }
        for (i=0L;i<ordermaxnumterms;i++)
        {
            free(S[i]);
        }
        free(S);
        
        if (nulldimflag==1)
        {
            break;
        }
        //Repeat with higher order
    } 
    
    
    for (i=0L;i<(MAX_ODE_ORDER+1L);i++)
    {
        for (j=0L;j<NUM_COEFFS;j++)
        {
            mpz_clear(D[i][j]);
        }
    }
    for (i=0L;i<(MAX_ODE_ORDER+1L);i++)
    {
        free(D[i]);
    }
    free(D);
    
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
        fprintf(fouteqs,"#ODE candidate for %s\n",finname);
        
        //for (n=0L;n<nulldim;n++)
        n=bestnulldim;
        {
            firstorder=0L;
            for (i=0L;i<numterms;i++)
            {
                mpz_set_ui(temp,0L);
                for (j=0L;j<MAX_POLY_ORDER+1L;j++) 
                {
                    mpz_abs(temp2,N[(i+j*numterms)*nulldim+n]);
                    mpz_add(temp,temp,temp2);
                }
                if (mpz_cmp_ui(temp,0L)>0L)
                {
                    firstorder++;
                    if (firstorder>1L)
                    {
                        fprintf(fouteqs,"+");
                        printf("+");
                    }
                    fprintf(fouteqs,"(");
                    printf("(");
                    firstterm=0L;
                    for (k=0L;k<MAX_POLY_ORDER+1L;k++)
                    {
                        if (mpz_cmp_ui(N[(i+k*numterms)*nulldim+n],0L)!=0L)
                        {
                            firstterm++;
                            if ((firstterm>1L)&&(mpz_cmp_ui(N[(i+k*numterms)*nulldim+n],0L)>0L))
                            {
                                fprintf(fouteqs,"+");
                                printf("+");
                                if (mpz_cmp_ui(N[(i+k*numterms)*nulldim+n],1L)!=0L)
                                {
                                    gmp_fprintf (fouteqs, "%Zd*", N[(i+k*numterms)*nulldim+n]);
                                    gmp_fprintf (stdout, "%Zd*", N[(i+k*numterms)*nulldim+n]);
                                }
                            }
                            else if ((firstterm>1L)&&(mpz_cmp_ui(N[(i+k*numterms)*nulldim+n],0L)<0L))
                            {
                                mpz_abs(temp,N[(i+k*numterms)*nulldim+n]);
                                if (mpz_cmp_ui(temp,1L)!=0L)
                                {
                                    gmp_fprintf(fouteqs, "%Zd*", N[(i+k*numterms)*nulldim+n]);
                                    gmp_fprintf(stdout, "%Zd*", N[(i+k*numterms)*nulldim+n]);
                                }
                                else
                                {
                                    fprintf(fouteqs,"-");
                                    printf("-");
                                }
                            }
                            else //first term
                            {
                                if (k>0)
                                {
                                    mpz_abs(temp,N[(i+k*numterms)*nulldim+n]);
                                    if (mpz_cmp_ui(temp,1L)!=0L)
                                    {
                                        gmp_fprintf(fouteqs, "%Zd*", N[(i+k*numterms)*nulldim+n]);
                                        gmp_fprintf(stdout, "%Zd*", N[(i+k*numterms)*nulldim+n]);
                                    }
                                    else if (mpz_cmp_ui(N[(i+k*numterms)*nulldim+n],0L)<0L)
                                    {
                                        fprintf(fouteqs,"-");
                                        printf("-");
                                    }
                                }
                                else
                                {
                                    gmp_fprintf(fouteqs, "%Zd", N[(i+k*numterms)*nulldim+n]);
                                    gmp_fprintf(stdout, "%Zd", N[(i+k*numterms)*nulldim+n]);
                                }
                            }
                            if (k>0L)
                            {
                                if (k==1)
                                {
                                    fprintf(fouteqs,"x");
                                    printf("x");
                                }
                                else
                                {
                                    fprintf(fouteqs,"x^%ld",k);
                                    printf("x^%ld",k);
                                }
                            }
                            if (MAX_FOUND_ORDER<k)
                            {
                                MAX_FOUND_ORDER=k;
                            }
                        }
                    }
                    fprintf(fouteqs,")");
                    printf(")");
                    for (j=0L;j<ODE_ORDER+1L;j++)
                    {
                        if (orderexp[i][j]!=0L)
                        {
                            if (j==0L)
                            {
                                if (orderexp[i][j]!=1L)
                                {
                                    fprintf(fouteqs,"*(y(x))^%ld",orderexp[i][j]);
                                    printf("*y^%ld",orderexp[i][j]);
                                }
                                else
                                {
                                    fprintf(fouteqs,"*y(x)");
                                    printf("*y");
                                }
                            }
                            else if (j==1L)
                            {
                                if (orderexp[i][j]!=1L)
                                {
                                    fprintf(fouteqs,"*(diff(y(x),x))^%ld",orderexp[i][j]);
                                    printf("*(Dx)^%ld",orderexp[i][j]);
                                }
                                else
                                {
                                    fprintf(fouteqs,"*diff(y(x),x)");
                                    printf("*Dx");
                                }
                            }
                            else
                            {
                                if (orderexp[i][j]!=1L)
                                {
                                    fprintf(fouteqs,"*(diff(y(x),x$%ld))^%ld",j,orderexp[i][j]);
                                    printf("*(Dx^%ld)^%ld",j,orderexp[i][j]);
                                }
                                else
                                {
                                    fprintf(fouteqs,"*diff(y(x),x$%ld)",j);
                                    printf("*Dx^%ld",j);
                                }
                            }
                        }
                    }
                }
            }
            fprintf(fouteqs,":\n");
            printf("=0\n");
            printf("Confidence level: %02ld%%\n",(long) floor((double) 100L-100L*nonzeroterms/(NUM_COEFFS-ODE_ORDER)));
        }
        
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
    
    for (i=0L;i<MAX_COEFFS;i++)
    {
        mpz_clear(temparray[i]);
    }
    mpz_clears(temp,temp2,coeff,NULL);
    
    time(&end);
    printf("\nEllapsed time %.fs\n",difftime(end,start));
    exit(EXIT_SUCCESS);
}
