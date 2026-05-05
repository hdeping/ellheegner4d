#29 ellheegner4d
## Abstract

This paper introduces a computational tool designed for computing
the generators of rank 1 elliptic curves. Based on a slightly modified
version of Pari/GP, the tool adopts a combined method of Heegner points
and 4-descent, which can effectively accelerate the computation of
generators. The above computational process is derived from Watkins'
work in 2006, and the computation of 4-descent makes use of the online
Magma computational tool. In particular, it should be noted that we
have used this tool to compute the generator of an elliptic curve
with a canonical height of 47349. The code of this paper has been
open-sourced on the GitHub platform.

Please find the complete article in the file computing47349.pdf.
## Demo of Computing the Generator

```
\\ height 47349
E = ellinit([0, -1, 1, -636884234186034, -6186412091784687672887]);
X = [x1,x2,x3,x4];
s1 = 18*x1^2 + 14*x1*x2 + 7*x1*x3 + 12*x1*x4 - 18*x2^2 - 40*x2*x3 - 20*x2*x4 + 28*x3^2 - 21*x3*x4 + 35*x4^2;
s2 = 17*x1^2 - 23*x1*x2 + 10*x1*x3 + 56*x1*x4 + 67*x2^2 + 14*x2*x3 +107*x2*x4 + 40*x3^2 - 35*x3*x4 - 34*x4^2;

M1 = get_mat_4(s1,X);
M2 = get_mat_4(s2,X);
height = ellL1(E,1)/ellbsd(E);
print("height = ",height);
\\ computing the point on the 4-cover
P = ellheegner_4descent(E,M1,M2,height);
\\ computing the generator
P1 = get_minimalmodel_point(M1,M2,P);

```

## Test Functions

The `test.gp` file contains a test function that demonstrates the usage of the ellheegner4d method to calculate generators for different elliptic curves. This test function showcases the application of the Heegner point and 4-descent method on various elliptic curves with different heights and coefficients.

The test function typically:
- Initializes an elliptic curve with specific coefficients
- Defines quadratic forms `s1` and `s2` to construct matrices `M1` and `M2`
- Computes the L-function value to determine the expected height
- Uses the `ellheegner_4descent()` function to compute the generator
- Converts the result to the minimal model using `get_minimalmodel_point()`

