Implementation:

	This is an implementation of paper "Memory efficient scheduling of Strassen-Winograd’s matrix multiplication algorithm"
	Refer to the paper's "Table 3: IP schedule for operation C <- A×B in place" for detailed steps.
	The paper could be downloaded from http://lig-membres.imag.fr/pernet/Publications/fp05-dumas.pdf

Compiler:

	The tested compiler is "Apple LLVM version 7.3.0 (clang-703.0.31)"

Code Structure:

	verify.c  : This is used to verify the implementation of "strassen.c".
	strassen.c : The implementation of the trassen-Winograd’s matrix multiplication algorithm.
	strassen.h : Header file of strassen.c.

Verification:

	Use the traditional sequential matrix multiplication to generate correct result for verification.

Efficiency comparison:

	Use the traditional sequential matrix multiplication to compare calculation time spent.


Usage :

	./verify opt1 opt2

	opt1 : The input square matrix's size N

	opt2 : Whether we want to verify and compare the calculation result with sequential matrix multiplication. 1:Yes  0:No

	Example : ./verify 2048 1


Note :

	Using ICC compiler might get better performance than clang. 
	[http://llvm.org/docs/Vectorizers.html]
