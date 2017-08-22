
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gmp.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"


//Compile: gcc -Wall serintode_flint.c -o serintode_flint.o -lflint -lgmp -lm -I /usr/local/include/flint
//There are issues with the results. It can give spurious results unless enough checks are used. This could maybe be caused by the lack of checking of the final result (?)
//Test case: file new_Chi3_w.ser, 20 checks, 400 MAX_COEFFS

int main()
{
    time_t start,end;
    char *finname = "tests/catalan.txt"; /*File name of data*/
    long const NUM_CHECKS=2L; /*Should be greater than 0*/
    long const MIN_ODE_ORDER=1L; 
    long const MAX_COEFFS=300; /*Should be checked for very large sequences*/
    long const MAX_LINE_LENGTH=100000L; 
    //char foutsumname[64]; /*Output summary file name*/
    char fouteqsname[64]; /*Output equations file name*/
    long NUM_COEFFS=0L;
    long MAX_ODE_ORDER=0L; 
    long ODE_ORDER=0L;
    long MAX_POLY_ORDER=0L;
    long MAX_FOUND_ORDER=0L;
    long COLUMNS=0L, ROWS=0L;
    long i,j,k,n,nonzeroterms,ordersused,termsused, MIN_MAX_FOUND_POLY_ORDER,MAX_FOUND_POLY_ORDER,MAX_FOUND_ODE_ORDER, firstorder,firstterm, mintermsused, bestnulldim, finalorders;
    long nulldim=0L;
    char input_string[MAX_LINE_LENGTH+1L];
    fmpz_t temp, temp2, coeff;
    fmpz_mat_t S, M, N;
    FILE *fin=NULL, *fouteqs=NULL; //, *foutsum=NULL
    char *fgcheck, nulldimflag=0;
    
    setvbuf(stdout,NULL,_IONBF,0);
    
    time(&start);
    
    
    fin = fopen(finname,"r");
    if (fin==NULL)
    {
        printf("\nError: Could not open input file %s. %s\n",finname,strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    fmpz_init(temp);
    fmpz_init(temp2);
    fmpz_init(coeff);
    fmpz_mat_init(S,1,MAX_COEFFS);  
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
                fmpz_mat_clear(S);
                fclose(fin);
                exit(EXIT_FAILURE);
            }
        }   
        else
        {
            fmpz_set_str(temp,input_string,10);
            if (nonzeroterms==0L) //Remove initial zeros
            {
                if (fmpz_cmp_ui(temp,0)!=0)
                {
                    nonzeroterms++;
                    fmpz_set(fmpz_mat_entry(S,0,NUM_COEFFS),temp);
                    NUM_COEFFS++;
                }
            }
            else
            {
                fmpz_set(fmpz_mat_entry(S,0,NUM_COEFFS),temp);
                NUM_COEFFS++;
            }
        }
    }
    fclose(fin); 
    
    /*
    printf("Input S Matrix:\n");
    for (i=0L;i<NUM_COEFFS;i++)
    {
        fmpz_print(fmpz_mat_entry(S,0,i)); printf(" ");
    }
    printf("\n\n");
    */
    
    fmpz_mat_init(M,MAX_COEFFS,MAX_COEFFS);
    printf("File %s:\n",finname);
    MAX_ODE_ORDER=floor((double) (NUM_COEFFS-NUM_CHECKS)/2L)-1L;
    for (n=MIN_ODE_ORDER;n<MAX_ODE_ORDER+1L;n++)
    {
        ODE_ORDER=n; //printf("Starting ODE order %ld\n",ODE_ORDER);
        MAX_POLY_ORDER=floor((double) (NUM_COEFFS-NUM_CHECKS-ODE_ORDER)/(ODE_ORDER+1L))-1L;
        if (MAX_POLY_ORDER == 0L)
        {
            break;
        }
        COLUMNS=(ODE_ORDER+1L)*(MAX_POLY_ORDER+1L);
        ROWS=COLUMNS+NUM_CHECKS; 
        printf("Checking for ODE order %ld with polynom coeffs of order %ld and %ld checks\n",ODE_ORDER,MAX_POLY_ORDER,NUM_CHECKS);
        fmpz_mat_init(M,ROWS,COLUMNS);
        
        //Null vector has the form: 0th order poly coeffs of first term up to ODE order term, order 1 poly coeffs of all terms, order 2 poly coeffs of all terms, etc
        //Matrix rows correspond to each successive order that needs to be annihilated
        //Matrix columns are the result of the derivatives and poly coeffs acting on input series to give a particular final order
        for (i=0L;i<MAX_POLY_ORDER+1L;i++)
        {
            for (j=0;j<ODE_ORDER+1L;j++)
            {
                for (k=i;k<ROWS;k++)
                {
                    fmpz_fac_ui(temp,j+k-i); //fmpz_print(temp); printf(" ");
                    fmpz_fac_ui(temp2,k-i); //fmpz_print(temp2); printf(" ");
                    fmpz_divexact(coeff,temp,temp2);  //coeff = (j+k-i)!/(k-i)!
                    fmpz_mul(temp,coeff,fmpz_mat_entry(S,0,j+k-i)); 
                    fmpz_set(fmpz_mat_entry(M,k,i*(ODE_ORDER+1L)+j),temp);
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
                    fmpz_print(fmpz_mat_entry(M,k,i*(ODE_ORDER+1L)+j)); printf(" ");
                }printf("\n");
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
                fmpz_print(fmpz_mat_entry(M,i,j)); printf(" ");
            }
            printf("\n");
        }
        */
        
        fmpz_mat_init(N,COLUMNS,COLUMNS);
        
        //nulldim = nullspaceMP(ROWS,COLUMNS,M,&N);
        nulldim = fmpz_mat_nullspace(N,M);
        
        fmpz_mat_clear(M);
        
        /*
        printf("Output N matrix\n");
        for (i=0L;i<COLUMNS*nulldim;i++)
        {
            gmp_fprintf (stdout, "  %Zd", N[i]);
        }
        printf("\n\n");
        */
        
        fmpz_set_ui(temp,0L);
        if (nulldim>0L)
        {
            //Checking whether kernel solution is spurious
            for (i = 0L; i < COLUMNS; i++)
            {
                for (j = 0L; j < nulldim; j++)
                {
                    fmpz_abs(temp2,fmpz_mat_entry(N,i,j));
                    fmpz_add(temp,temp,temp2);
                }
            }
            if (fmpz_get_ui(temp)!=nulldim)
            {
                nulldimflag=1;
            }
        }
        
        if (nulldimflag==1)
        {
            //Check that the number of orders is not 1
            mintermsused=0;
            MIN_MAX_FOUND_POLY_ORDER=MAX_POLY_ORDER+1L;
            bestnulldim=0L;
            finalorders=0;
            for (k=0;k<nulldim;k++)
            {
                ordersused=0L;
                termsused=0L;
                nonzeroterms=0L;
                MAX_FOUND_POLY_ORDER=0L;
                MAX_FOUND_ODE_ORDER=0L;
                for (i=0L;i<ODE_ORDER+1L;i++)
                {
                    for (j=0L;j<MAX_POLY_ORDER+1L;j++)
                    {
                        if (fmpz_cmp_ui(fmpz_mat_entry(N,i+j*(ODE_ORDER+1L),k),0L)!=0L)
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
                        ordersused++;
                        MAX_FOUND_ODE_ORDER=i;
                    }
                }
                if (MAX_FOUND_ODE_ORDER==ODE_ORDER)
                {
                    if (MIN_MAX_FOUND_POLY_ORDER>MAX_FOUND_POLY_ORDER)
                    {
                        MIN_MAX_FOUND_POLY_ORDER=MAX_FOUND_POLY_ORDER;
                        mintermsused=termsused;
                        finalorders=ordersused;
                        bestnulldim = k;
                    }
                    else if (MIN_MAX_FOUND_POLY_ORDER==MAX_FOUND_POLY_ORDER)
                    {
                        if (mintermsused>termsused)
                        {
                            mintermsused=termsused;
                            finalorders=ordersused;
                            bestnulldim = k;
                        }
                    }
                }
            }
            if (finalorders<2L)
            {
                //printf("Spurious equation with %ld order term.\n",finalorders);
                nulldimflag=0;
                fmpz_mat_clear(N);
            }
            else if (mintermsused<finalorders+1L)
            {
                //printf("Polynomial coefficients only have one term each.\n");
                nulldimflag=0;
                fmpz_mat_clear(N);
            }
            else if (MAX_ODE_ORDER<ODE_ORDER)
            {
                //printf("Spurious solution came from taking a lot of derivatives, then found lower order ODE.\n");
                nulldimflag=0;
                fmpz_mat_clear(N);
            }
            else
            {
                break;
            }
        }
        else
        {
            fmpz_mat_clear(N);
        }
        
    }
    
    fmpz_mat_clear(S);
    
    
    
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
        printf("\n***********************\n");
        printf("***Found a solution!***\n");
        printf("*Confidence level: %02ld%%*\n",(long) floor((double) 100L-100L*(ODE_ORDER+1L)*(MIN_MAX_FOUND_POLY_ORDER+1L)/(NUM_COEFFS-ODE_ORDER)));
        printf("***********************\n\n");
        sprintf(fouteqsname,"%s_solution_%ld-checks.txt",finname,NUM_CHECKS);
        fouteqs = fopen(fouteqsname,"w");
        if (fouteqs==NULL)
        {
            printf("\nError: Could not open equations output file %s. %s\n",fouteqsname,strerror(errno));
            fclose(fin);
            fmpz_mat_clear(N);
            exit(EXIT_FAILURE);
        }
        setvbuf(fouteqs,NULL,_IOLBF,32);
        
        
        i=0L;
        while (fmpz_cmp_ui(fmpz_mat_entry(N,i,bestnulldim),0L)!=0L)
        {
            i++;
        }
        fmpz_set(coeff,fmpz_mat_entry(N,i,bestnulldim));
        for (j=i+1L;j<COLUMNS;j++)
        {
            if (fmpz_cmp_ui(fmpz_mat_entry(N,j,bestnulldim),0L)!=0L)
            {
                fmpz_gcd(coeff,coeff,fmpz_mat_entry(N,j,bestnulldim));
            }
        } //printf("GCD = "); fmpz_print(coeff); printf("\n");
        firstorder=0L;
        fprintf(fouteqs,"ODE%s := ",finname);
        for (i=0L;i<ODE_ORDER+1L;i++)
        {
            fmpz_set_ui(temp,0L);
            for (j=0L;j<MAX_POLY_ORDER+1L;j++)
            {
                fmpz_abs(temp2,fmpz_mat_entry(N,(i+j*(ODE_ORDER+1L)),bestnulldim));
                fmpz_add(temp,temp,temp2);
            }
            if (fmpz_cmp_ui(temp,0L)>0L)
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
                    if (fmpz_cmp_ui(fmpz_mat_entry(N,(i+k*(ODE_ORDER+1L)),bestnulldim),0L)!=0L)
                    {
                        firstterm++;
                        if ((firstterm>1L)&&(fmpz_cmp_ui(fmpz_mat_entry(N,(i+k*(ODE_ORDER+1L)),bestnulldim),0L)>0L))
                        {
                            fprintf(fouteqs,"+");
                            printf("+");
                            fmpz_divexact(temp,fmpz_mat_entry(N,(i+k*(ODE_ORDER+1L)),bestnulldim),coeff);
                            if (fmpz_cmp_ui(temp,1L)!=0L)
                            {
                                fmpz_fprint(fouteqs,temp); fprintf(fouteqs,"*");
                                fmpz_fprint(stdout, temp); fprintf(stdout,"*");
                            }
                        }
                        else if ((firstterm>1L)&&(fmpz_cmp_ui(fmpz_mat_entry(N,(i+k*(ODE_ORDER+1L)),bestnulldim),0L)<0L))
                        {
                            fmpz_divexact(temp,fmpz_mat_entry(N,(i+k*(ODE_ORDER+1L)),bestnulldim),coeff);
                            fmpz_abs(temp2,temp);
                            if (fmpz_cmp_ui(temp2,1L)!=0L)
                            {
                                fmpz_fprint(fouteqs,temp); fprintf(fouteqs,"*");
                                fmpz_fprint(stdout, temp); fprintf(stdout,"*");
                            }
                            else
                            {
                                fprintf(fouteqs,"-");
                                printf("-");
                            }
                        }
                        else
                        {
                            fmpz_divexact(temp,fmpz_mat_entry(N,(i+k*(ODE_ORDER+1L)),bestnulldim),coeff);
                            if (k>0)
                            {
                                fmpz_abs(temp2,temp);
                                if (fmpz_cmp_ui(temp2,1L)!=0L)
                                {
                                    fmpz_fprint(fouteqs,temp); fprintf(fouteqs,"*");
                                    fmpz_fprint(stdout, temp); fprintf(stdout,"*");
                                }
                                else if (fmpz_cmp_ui(temp,0L)<0L)
                                {
                                    fprintf(fouteqs,"-");
                                    printf("-");
                                }
                            }
                            else
                            {
                                fmpz_fprint(fouteqs,temp);
                                fmpz_fprint(stdout, temp);
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
                if (i>0L)
                {
                    if (i==1)
                    {
                        fprintf(fouteqs,"*diff(y(x),x)");
                        printf("*Dx");
                    }
                    else
                    {
                        fprintf(fouteqs,"*diff(y(x),x$%ld)",i);
                        printf("*Dx^%ld",i);
                    }
                }
                else
                {
                    fprintf(fouteqs,"*y(x)");
                    printf("*y");
                }
            }
        }
        fprintf(fouteqs,":\n");
        printf("=0\n");
        fmpz_mat_clear(N);
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

    fmpz_clear(temp);
    fmpz_clear(temp2);
    fmpz_clear(coeff);
    return 0;
}
