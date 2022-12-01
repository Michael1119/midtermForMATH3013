# Midterm for MATH3013 - A first order ODE solver

## Features
* A Runge-Kutta class supporting up to 4th order Runge-Kutta method
* A first order ODE class including Lorenz system and Newton's cooling law

## Methods
lorenz(string filename) / cooling(string filename) : create a Lorenz / Newton's cooling law system
* .derivative : return the derivative of the system
* .get_t0 : return start time
* .get_tn : return end time
* .get_h : return step size
* .get_n : return number of steps
* .get_x : return initial condition
* .get_filename : return output filename

RungeKutta(int order) : create a RK method with a specific order of local truncation error (up to 4)
* .solve(lorenz / cooling) : solve the first order ODE system and save the results to the output file