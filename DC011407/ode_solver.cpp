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
    k = 0.03; // heat transfer coefficient
    Ts = 20.; // surrounding temperature
}

void cooling::setupDomainAndPartition(){
    t0 = 0.; // start time
    tn = 40.; // end time
    h = 0.01; // step size
    n = tn / h; // number of steps
}

void cooling::setupInitialCondition(){
    x = {40.}; // initial temperature
}