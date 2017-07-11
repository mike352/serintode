serintode is a C program for searching for ODEs which annihilate an integer series up to a certain order using the arbitrary precision integer library IML which is based on GMP.












OLD VERSION THAT RELIES ON IML

serintode_iml uses IML, the Integer Matrix Library, https://cs.uwaterloo.ca/~astorjoh/iml.html
IML uses GMP, but not fully. The functions we make use of from IML are kernelMP or nullspaceMP. These functions reduces the input matrix mod p, where p is a random prime. From nullspace.c, we see that these functions both have the following line:
p = RandPrime(15, 19);
The function RandPrime is defined in basisop.c, where the statement above means that the random prime is found between 2^15 and 2^19-1. 

If any input matrix element is 0 mod p, then it assigns it the value p according to this line in nullspace.c:
DA[i] = (double) ((temp = (A[i] % ((long) p))) >= 0 ? temp : ((long) p) + temp);
This doesn’t give the correct result if the integer is large enough and has multiple factors of p. What is more, we can see in the line above that the result is converted to a double (!) because IML it uses CBLAS in its RowEchelonTransform function. 

A fix isn't simple. The RowEchelonTransform function has the precondition:
ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
where n is the number of rows. This seems to all stem from the fact that it reduces everything mod p and then uses a double type. So artificially increasing the range of the random primes cannot help, the RowEchelonTransform function would have to be rewritten to not rely on CBLAS.

So as long as none of series coefficients have multiple factors of the random prime chosen, the results can be trusted. Otherwise, not. If the coefficients are large enough, the serintode_iml function should be run multiple times and any result should be checked explicitly using the original series to check for annihilation. This could be implemented at a later time.

Therefore, I'm switching over to another fully arbitrary precision integer matrix library.
