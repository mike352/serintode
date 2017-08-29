#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "gmp.h"

//Only used to create a data directory or check that one already exists
 #if (defined(__unix__) || defined(__APPLE__))
	#include <sys/stat.h> 
	#define OS 1 
	int _mkdir(char*); //prototyping to make warning go away when compiling with gcc -Wall
#elif (defined(_WIN16) || defined(_WIN32)|| defined(_WIN64))
	#include <direct.h> 
	#define OS 0
#else //Otherwise ask
	#define OS -1 
#endif

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


int main(int argc, char* argv[])
{
    long const MIN_ODE_ORDER=1L;
    long const MIN_DEPTH=1L;
    long MAX_ODE_ORDER,NUM_COEFFS,ODE_ORDER,MAX_DEPTH_POSS;
    long i,j,n,p, *s, **orderexp, arrindex,numterms,maxnumterms;
    char fname[64], *dirname="lookuptables";
    FILE *fout=NULL;
    mpz_t temp,temp2;
    int dcheck;
    
    setvbuf(stdout,NULL,_IONBF,0);
    
    if (argc!=3)
    {
        printf("\nWrong usage: Pass two arguments, max order and number of series coefficients. Example: ./makelookup 100 10000\n");
        exit(EXIT_FAILURE);
    }
    MAX_ODE_ORDER = atol(argv[1]);
    NUM_COEFFS = atol(argv[2]);
    maxnumterms = NUM_COEFFS;
    
   
    //Create directory if necessary
    if (OS) //UNIX
    {
        dcheck = mkdir(dirname,S_IRWXU | S_IRWXG | S_IRWXO);
        if ((dcheck==-1)&(errno!=EEXIST))
        {
            printf("\nERROR: Could not create output directory. %s\n",strerror(errno));
            exit(EXIT_FAILURE);
        }
    }
    else if (OS==0) //Windows
    {
        dcheck = _mkdir(dirname);
        if ((dcheck==-1)&(errno!=EEXIST))
        {
            printf("\nERROR: Could not create output directory. %s\n",strerror(errno));
            exit(EXIT_FAILURE);
        }
    }
    else if (OS==-1) //Ask instead
    {
        printf("\nPlease ensure a directory called ""lookuptables"" exists already.\n");
    }
   
    
    mpz_inits(temp,temp2,NULL);
    for (n=MIN_ODE_ORDER;n<=MAX_ODE_ORDER;n++)
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
        }
        for (p=MIN_DEPTH;p<=MAX_DEPTH_POSS;p++)
        {
            mpz_fac_ui(temp,ODE_ORDER+1L+p);
            mpz_fac_ui(temp2,ODE_ORDER+1L);
            mpz_divexact(temp,temp,temp2);
            mpz_fac_ui(temp2,p);
            mpz_divexact(temp,temp,temp2);
            numterms=mpz_get_ui(temp)-1L; //Equals (n+p)!/(n!p!) Need one less than this
            
            s = (long *) calloc((ODE_ORDER+1),sizeof(long));
            orderexp = (long **) malloc((numterms+1L)*sizeof(long *));
            if ((orderexp==NULL) || (s==NULL))
            {
                fprintf(stderr, "No memory left for allocating orderexp matrix. %s",strerror(errno));
                exit(EXIT_FAILURE);
            }
            for (i=0;i<(numterms+1L);i++)
            {
                orderexp[i] = (long *) calloc((ODE_ORDER+1),sizeof(long));
                if(orderexp[i]==NULL)
                {
                    fprintf(stderr, "No memory left for allocating orderexp array. %s",strerror(errno));
                    for (j=0L;j<i;j++)
                    {
                        free(orderexp[j]);
                    }
                    free(s);
                    exit(EXIT_FAILURE);
                }
            }
            arrindex=0L;
            
            combs(orderexp,s,ODE_ORDER+1L,p,0,&arrindex);
            
            sprintf(fname,"%s/o%ldd%ld.txt",dirname,n,p);
            fout = fopen(fname,"w");
            if (fout==NULL)
            {
                printf("Couldn't open the output file %s. %s\n",fname,strerror(errno));
                exit(EXIT_FAILURE);
            }
            for (i=0L;i<numterms;i++)
            {
                for (j=0L;j<=ODE_ORDER;j++)
                {
                    fprintf(fout,"%ld ",orderexp[i][j]);
                }
                fprintf(fout,"\n");
            }
            fclose(fout);
                
            
            for (i=0L;i<(numterms+1L);i++)
            {
                free(orderexp[i]);
            }
            free(orderexp);
            free(s);
            
        } //End depth p loop
    } //End order n lopo
    
    mpz_clears(temp,temp2,NULL);
    exit(EXIT_SUCCESS);
}