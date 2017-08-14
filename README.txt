serintode is a C program for searching for ODEs which annihilate an integer series up to a certain order.

There are two programs based on two different libraries, IML and flint. The IML version has linear and algebraic ODE search versions. At the moment the flint version only has a linear ODE version. 

Both versions rely on the construction of a matrix nullspace, the solution being one of the null vectors (all are solutions). Due to the way each library constructs the nullspace, the output from each program can be different. Both outputs give the corred ODE. The flint solutions may have smaller polynomial coefficient orders, although I haven't tested the differences extensively enough to say that's always or even usually the case.

The program automatically determines the number of coefficients in a file and searches for ODEs of increasing order according to the number of coefficients and the number of checks required. It outputs the result both to the screen as well as a file created if a solution is found, using the input file name as the start of the file.

Input files should have one coefficient per line.

In the current version the C files should be modified for changing the input file name, the number of checks, and the minimum ODE order. 

If there are more than 10,000 coefficients, MAX_COEFFS should be changed accordingly.

If the coefficients are larger than 100,000 digits, the MAX_LINE_LENGTH should be changed accordingly.

The algebraic ODE version serintode_iml_nonlin.c searches for algebraic ODEs up to an order determined by the number of coefficients and up to n-products with the "depth" n determined also by the number of coefficients. An order 2 depth 3 ODE would have terms:
p1*y + p2*Dx + p3*Dx^2 + p4*y^3 + p5*y^2*Dx + p6*y^2*Dx^2 + p7*y*(Dx)^2 + p8*y*Dx*Dx^2 + p9*y*(Dx^2)^2 + p10*(Dx)^3 + p11*(Dx)^2*Dx^2 + p12*Dx*(Dx^2)^2 + p13*(Dx^2)^3
where p_k are polynomial coefficients.


VERSIONS THAT RELY ON IML
serintode_iml uses IML, the Integer Matrix Library, https://cs.uwaterloo.ca/~astorjoh/iml.html, which relies on GMP, CBLAS, and possibly ATLAS.
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

The IML version is faster than the flint version.

Compilation can be done via:
gcc -Wall serintode_iml.c -o serintode_iml.o -liml -lcblas -lgmp -lm
gcc -Wall serintode_iml_nonlin.c -o serintode_iml_nonlin.o -liml -lcblas -lgmp -lm


VERSION THAT RELIES ON flint
serintode_flint uses the flint library, http://www.flintlib.org/, which uses fully arbitrary precision integers via GMP for all of its calculations, without using intermediary mod p calculations. In testing, it is slower than IML. For simply determining whether to construct the nullspace, that is, whether COLUMNS-RANK>0, it is only roughly 5% slower. But for the actual construction of the nullspace, it can be 3x slower for a matrix of size 400 or so. IML checks whether the matrix annihilates the nullspace vectors, but flint doesn't appear to do that. I plan to add that in myself. That may be the reason why flint sometimes gives spurious results unless enough checks are given. More testing is needed to determine that question. 

Compilation can be done via:
gcc -Wall serintode_flint.c -o serintode_flint.o -lflint -lgmp -lm -I /usr/local/include/flint