# Discontinuous Galerkin
 This repository contains the _operators_ which form an integral part of the **Discontinuous Galerkin methods**. The focus is to make it usable for fluid dynamics.


## Functions
 - Essentials
    - `zeros(double *A,unsigned m, unsigned n)`
        - **Dependencies**: _None_
        - **Input**:
            -  `double *A` : _The pointer to the first element of the matrix_.
            -  `unsigned m` : _The number of rows of the matrix_.
            -  `unsigned n` : _The number of columns of the matrix_.  
        - **Brief**: This is used to create a matrix full of zeros. Can be used as a basic initialization function. _Caution_: The matrix should already be allocated memory, that is the assumption of this subroutine.
        - **Usage**: `zeros(&A[0][0],6,3)`, where `A` is a statically declared 2-D array with dimensions `6x3`.
    - `ones(double *A,unsigned m, unsigned n)`
        - **Dependencies**: _None_
        - **Input**:
            -  `double *A` : _The pointer to the first element of the matrix_.
            -  `unsigned m` : _The number of rows of the matrix_.
            -  `unsigned n` : _The number of columns of the matrix_.  
        - **Brief**: This is used to create a matrix full of ones. Can be used as a basic initialization function. _Caution_: The matrix should already be allocated memory, that is the assumption of this subroutine.
        - **Usage**: `ones(&A[0][0],6,3)`, where `A` is a statically declared 2-D array with dimensions `6x3`.
    - `display(double *A,unsigned m, unsigned n)`
        - **Dependencies**: _None_
        - **Input**:
            -  `double *A` : _The pointer to the first element of the matrix_.
            -  `unsigned m` : _The number of rows of the matrix_.
            -  `unsigned n` : _The number of columns of the matrix_.  
        - **Brief**: This is used a display a matrix of dimensions `mxn`. The default function would display 2 floating point numbers but this can be easily _tweaked_.
        - **Usage**: `display(&A[0][0],6,3)`, where `A` is a statically declared 2-D array with dimensions `6x3`.


## Contact

If you have any questions about the usage or contribution, drop me a mail at `kaushikcfd@gmail.com`.
