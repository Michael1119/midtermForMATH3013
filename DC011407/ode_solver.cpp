#include "ode_solver.h"

void lorenz::setupParameter(){
    sigma = 10.;
    R = 28.;
    b = 10./3;
}

void lorenz::setupDomainAndPartition(){
    t0 = 0.; // start time
    tn = 40.; // end time
    h = 0.01; // step size
    n = tn / h; // number of steps
}

void lorenz::setupInitialCondition(){
    x = {1., 1., 1.}; // initial position
}

void cooling::setupParameter(){
    k = 0.1; // heat transfer coefficient
    Ts = 293.; // surrounding temperature (unit : K)
}

void cooling::setupDomainAndPartition(){
    t0 = 0.; // start time
    tn = 10.; // end time
    h = 2; // step size
    n = tn / h; // number of steps
}

void cooling::setupInitialCondition(){
    x = {313.}; // initial temperature (unit : K)
}