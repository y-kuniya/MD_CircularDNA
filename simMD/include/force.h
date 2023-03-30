#ifndef __FORCE_H__
#define __FORCE_H__

#include"variables.h"
#include"systemparam.h"

class Force{
    private:
    int PN;
    int SIZE;
    vector<double> cosBeta;
    vector<double> Beta_over_sinBeta;
    void clear(void);
    void calc_excluedVolume(Variables *vars);
    void calc_stretching(Variables *vars);
    void calc_Beta(Variables *vars);
    void calc_bending(Variables *vars);

    public:
    vector<double> F;
    void init(const int PN_);
    void calc_force(Variables *vars);
};

#endif