# Discontinuous Galerkin
 This repository contains the _operators_ which form an integral part of the **Discontinuous Galerkin methods**. The focus is to make it usable for fluid dynamics.


## Functions
 - Essentials
    - `zeros(double **A,unsigned m, unsigned n)`
        - **Dependencies**: _None_
        - **Input**:
            -  `double **A` : _The address of the matrix_.
            -  `unsigned m` : _The number of rows of the matrix_.
            -  `unsigned n` : _The number of columns of the matrix_.  
        - **Brief**: This is used to create a matrix full of zeros. Can be used as a basic initialization function. _Caution_: The matrix should already be allocated memory, that is the assumption of this subroutine.
    - `ones(double **A,unsigned m, unsigned n)`
        - **Dependencies**: _None_
        - **Input**:
            -  `double **A` : _The address of the matrix_.
            -  `unsigned m` : _The number of rows of the matrix_.
            -  `unsigned n` : _The number of columns of the matrix_.  
        - **Brief**: This is used to create a matrix full of ones. Can be used as a basic initialization function. _Caution_: The matrix should already be allocated memory, that is the assumption of this subroutine.


## Contact

If you have any questions about the usage or contribution, drop me a mail at `kaushikcfd@gmail.com`.
