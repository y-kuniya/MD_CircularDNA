#ifndef __OBSERVER_H__
#define __OBSERVER_H__

#include"variables.h"

class Observer{
    private:
    double c0[3];
    public:
    void set_c0(Variables *vars);
    double squared_distance(Variables *vars);
};

#endif