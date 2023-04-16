#ifndef __MD_DRY_H__
#define __MD_DRY_H__

#include "variables.h"
#include "observer.h"
#include "systemparam.h"
#include "tools.h"
#include "force.h"


class MDdry{
    private:
    Variables *vars;
    Variables *vars_pred;
    Observer  *obs ;
    Force     *force;
    vector<double> S;
    vector<double> DT;
    vector<double> u_sub;
    int PN;
    int SIZE;
    void set_PN(int PN_);
    void makeconf(void);
    void update_pred(vector<double> &R);
    void update_corr(vector<double> &R);
    void update_frame_pred(vector<double> &R);
    void update_frame_corr(vector<double> &R);

    public:
    MDdry(void);
    ~MDdry(void);
    void make_trajectory(void);
};




#endif