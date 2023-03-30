#include "../include/observer.h"

void 
Observer::set_c0(Variables *vars){
    vars->calc_center();
    for(int i=0;i<3;i++){
        c0[i] = vars->c[i];
    }
}

double
Observer::squared_distance(Variables *vars){
    double sd = 0.0;
    vars->calc_center();
    for(int i=0;i<3;i++){
        double delta = vars->c[i] - c0[i];
        sd += delta*delta;
    }
    return sd;   
}
