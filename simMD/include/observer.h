#ifndef __OBSERVER_H__
#define __OBSERVER_H__

#include"variables.h"
#include"force.h"
#include <fstream>

class Observer{
    private:
    double c0[3];
    ofstream ofs;

    public:
    Observer(void);
    ~Observer(void);
    void set_c0(Variables *vars);
    double squared_distance(Variables *vars);
    double writhing_number(Variables *vars);
    double total_twist(Force *force);
    void export_Wr_Tw(int step,Variables *vars,Force *force);
};

#endif