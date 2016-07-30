# Discontinuous Galerkin

## Introduction
 This repository contains the _operators_ which form an integral part of the **Discontinuous Galerkin methods**. The focus is to make it usable for fluid dynamics. The current code is able to solve the Shallow Water Hyperbolic set of equations in a two-dimensional space.

 The result of the above code for a Gaussian Initial Input, after visualizing in ParaView are as follows:

[![Solution using N = 4](http://img.youtube.com/vi/TnMSwGVz8CQ/0.jpg)](http://www.youtube.com/watch?v=TnMSwGVz8CQ)

## Advantages of the Solver
* **Uses BLAS and LAPACK librarires** . The problem had been broken down into various Linear Algebra problems like Matrix-Vector, Matrix-Matrix Multiplication and Gaussian Elimination. To mimimise the computational run-time, linear algebra sub-routines from BLAS and LAPACK were used.
* **Object Oriented** . The code was created so that it could be easily modififed for further development.
* **Ability to handle dicontinuous Inital Conditions**. The following simulation of Dam Break will illustrate it:

[![Solution using N = 4](http://img.youtube.com/vi/IZzkMfUm9OE/0.jpg)](http://www.youtube.com/watch?v=IZzkMfUm9OE)


## License

You can find all the details over here [LICENSE](https://github.com/kaushikcfd/Discontinuous-Galerkin/blob/master/LICENSE.md).

## Contact

If you have any questions about the usage or contribution, feel free to drop a mail at `kaushikcfd@gmail.com` or `sgopalak@iitb.ac.in`.
