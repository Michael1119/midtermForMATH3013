#include "ode_solver.h"

void lorenz::setupParameter(){
    sigma = 10.;
    R = 28.;
    b = 10./3;
}

void lorenz::setupDomainAndPartition(){
    t0 = 0.;
    tn = 40.;
    h = 0.01;
    n = tn / h;
}

void lorenz::setupInitialCondition(){
    x = {1., 1., 1.};
}

void cooling::setupParameter(){
    k = 0.03;
    Ts = 20.;
}

void cooling::setupDomainAndPartition(){
    t0 = 0.;
    tn = 40.;
    h = 0.01;
    n = tn / h;
}

void cooling::setupInitialCondition(){
    x = {40.};
}