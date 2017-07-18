serintode is a C program for searching for ODEs which annihilate an integer series up to a certain order using the arbitrary precision integer library IML which is based on CBLAS and GMP.

There are two programs, which will eventually be merged into one. The most current one is serintode_iml_auto.c, which automatically determines the number of coefficients in a file and searches for ODEs of increasing order according to the number of coefficients and the number of checks required. It outputs the result both to the screen as well as a file created if a solution is found, using the input file name as the start of the file.

serintode_iml_auto.c is recommended. The current serintode_iml.c outputs to the screen as follows.
If it finds a solution, it will output it as column vectors spanning the kernel of the input matrix. Each column vector has rows as follows:
P_{0,0}, P_{1,0}, ..., P_{o,0}, P_{0,1}, P_{1,1}, ..., P_{o,1}, ...  P_{0,n}, P_{1,n}, ..., P_{o,n}
where the order of the ODE is o and the order of the polynomial coefficients is n. The P_{i,j} are the j polynomial coefficients of the operator d^i/dx^i

Input files should have one coefficient per line.

In the current version the C files should be modified for changing the input file name, the number of checks, and the minimum ODE order. 

If there are more than 10,000 coefficients, MAX_COEFFS should be changed accordingly.

If the coefficients are larger than 100,000 digits, the MAX_LINE_LENGTH should be changed accordingly.

The current version requires IML, which requires GMP, CBLAS and possibly ATLAS. After installing dependencies, compilation can be done via
gcc -Wall serintode_iml.c -o serintode_iml.o -liml -lcblas -lgmp -lm
gcc -Wall serintode_iml_auto.c -o serintode_iml_auto.o -liml -lcblas -lgmp -lm
 


VERSION THAT RELIES ON IML

serintode_iml uses IML, the Integer Matrix Library, https://cs.uwaterloo.ca/~astorjoh/iml.html
IML uses GMP, but not fully. The functions we make use of from IML are kernelMP or nullspaceMP. These functions reduces the input matrix mod p, where p is a random prime. From nullspace.c, we see that these functions both have the following line:
p = RandPrime(15, 19);
The function RandPrime is defined in basisop.c, where the statement above means that the random prime is found between 2^15 and 2^19-1. 

If any input matrix element is 0 mod p, then it assigns it the value p, performing the mod p procedure according to this line in nullspace.c:
DA[i] = (double) ((temp = (A[i] % ((long) p))) >= 0 ? temp : ((long) p) + temp);
We see in the line above that the result is converted to a double, which is why the random prime needs to be in the range 2^15 and 2^19-1, or rather, as stated in the RowEchelonTransform function, it has the precondition:
ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
where n is the number of rows.
This is because IML uses CBLAS in its RowEchelonTransform function, which requires double types.

If IML finds a solution, it checks it, so any solution returned by IML is correct. However, if the check fails, it will attempt a new random prime and repeat the process. In principle, this could go on for a while, but in practice the input matrix would have to be quite pathological. A particular prime p is bad iff the gcd of all minors are divisible by p. Since in this case IML chooses another random prime, the condition that the input matrix will totally fail the algorithm is if the gcd of all minors is divisible by all primes in the range 2^15 and 2^19-1. In practice, if the gcd of all minors is divisible by a sufficient number of these primes, the algorithm could take a long while to find a suitable prime. This would still be rare. The theoretical algorithm is described here:
https://cs.uwaterloo.ca/~astorjoh/jscdense.pdf

I will be switching over to another fully arbitrary precision integer matrix library soon, flint, as a means of double-checking everything as well as to avoid any pathological cases.

