# PHSX815_HW7

This program integrates the function sin(pi\*x) + 1 on a closed interval, [0,1] with three separate methods, Monte Carlo integration (implemented as rejection sampling), Simpson's rule (n=2), and Gauss-Legendre quadrature (n=2). The user may input the number of MC steps to sample with.


-Nsteps					(integer) Sets the number of sampling steps in the Monte Carlo integration


The program prints out the analytic solution, the numerical solution, and their differences for each integration method. It also plots the Monte Carlo integral as a function of the number of samples, showing whether or not the integral converges to the analytic value and the estimates from the other two integration methods.

The program requires numpy and matplotlib.
