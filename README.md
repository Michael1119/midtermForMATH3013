# Midterm for MATH3013 - A first order ODE solver

## Features
* A first order ODE class including Lorenz system and Newton's cooling law
* A Runge-Kutta class supporting up to 4th order Runge-Kutta method

## Usage
Change the parameters in "ode_solver.cpp", then execute the following command in the terminal:

` g++ main.cpp ode_solver.cpp -o ode_solver `

This will generate an executable file. When executed, several .dat file will be generated.

The solution of the ODE system can be plotted by gnuplot.

Lorenz system :
```
set xlabel 'x1'
set ylabel 'x3'
plot 'lorenz-rk4.dat' using 2:4 with line
```
Newton's cooling law :
```
set xlabel 'Time (s)'
set ylabel 'Temperature (K)'
plot 'cooling-rk1.dat' with line,\
'cooling-rk2.dat' with line,\
'cooling-rk3.dat' with line,\
'cooling-rk4.dat' with line
```



## List of methods
lorenz(string filename) / cooling(string filename) : create a Lorenz / Newton's cooling law system
* .derivative(double time, valarray unknowns) : return the derivative of the system
* .get_t0() : return start time
* .get_tn() : return end time
* .get_h() : return step size
* .get_n() : return number of steps
* .get_x() : return initial condition
* .get_filename() : return output filename (default is "lorenz" / "cooling")

RungeKutta(int order) : create a RK method with a specific order of local truncation error (up to 4)
* .solve(lorenz / cooling) : solve the first order ODE system and save the results to the output file