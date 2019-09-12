README file for MATLAB code to accompany Vazquez, A. R. and Xu, H. (2019). "Construction 
of Two-Level Nonregular Designs of Strength Three With Large Run Sizes." Technometrics, 
Vol. 61, p.p. 341-353.

===CONTENTS=======================================

A. SCRIPT FILES
B. PRIMARY FUNCTIONS
C. UTILITY FUNCTIONS
D. DATA FILES

==================================================

A) SCRIPT FILES.
Use these to construct nonregular designs using the VNS algorithm.

Example_basic_prime.m       Shows how to construct a concatenated design from 
                            a regular design with a prime number of basic factors.
                 
Example_basic_no_prime.m    Shows how to construct a concatenated design from 
                            a regular design with a number of basic factors that
                            is NOT a prime.


B) PRIMARY FUNCTIONS.
These are the functions one would interact with when constructing and evaluating 
the designs. 

VNS.m                   Main function for the VNS algorithm for regular designs 
                        with a prime number of basic factors.
                        
VNS_noprime.m           Main function for the VNS algorithm for regular designs 
                        with a number of basic that is NOT a prime.
                                               
F4.m                    Computes the F4 vector and B4 value of a two-level design  
                        of strength 3.

rankX2.m                Computes the rank of the two-factor interaction matrix 
                        of a two-level strength-3 design.                        

C) UTILITY FUNCTIONS.
These are called by the primary functions for specific purposes.

vnsalgorithmfour.m      VNS algorithm for optimizing the F4 vector.

Concatenate.m           Concatenate several designs.
                        
Objective.m             Calculate F4 vector over selected 4-factor sets.
                        
FindSetConfounded.m     Find completely aliased words of length 4.
                        
D) DATA FILES.

regular_designs         Minimum aberration regular designs with two levels. The designs were
                        obtained from the FrF2 package in R. See Supplementary Section D for
                        details.

