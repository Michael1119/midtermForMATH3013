#include "ode_solver.h"

int main(){
    RungeKutta rk1(1);
    RungeKutta rk2(2);
    RungeKutta rk3(3);
    RungeKutta rk4(4);
    lorenz myLorenz;
    rk4.solve(myLorenz);
    cooling myCooling;
    rk1.solve(myCooling);
    rk2.solve(myCooling);
    rk3.solve(myCooling);
    rk4.solve(myCooling);
    return 0;
}