#ifndef __MD_H__
#define __MD_H__

#include "variables.h"
#include "observer.h"
#include "systemparam.h"
#include "tools.h"
#include "force.h"
#include <map>

class MD{
    private:
    Variables *vars;
    Variables *vars_pred;
    Observer  *obs ;
    Force     *force;
    SquareVec D_pred;
    SquareVec D_corr;
    SquareVec D_sqrt;
    vector<double> S;
    int PN;
    int SIZE;
    void set_PN(int PN_);
    void makeconf(void);
    void update_pred(vector<double> &R);
    void update_corr(vector<double> &R);
    public:
    MD(void);
    ~MD(void);
    void make_trajectory(void);
    void calc_sd(map<double,double> &m_sd,
                    const int trial,
                    const int STEP,
                    const int DISPOSE,
                    const int INTERVAL);
};

void run_m_sd(void);

#endif