#ifndef __VARIABLES_H__
#define __VARIABLES_H__

#include<vector>
using namespace std;

class Variables{
    public:
    vector<double> r;
    vector<double> f;
    vector<double> v;
    vector<double> u;
    vector<double> b;
    vector<double> inv_b;
    double num ;
    double c[3];

    Variables(void);
    void init(int num_) ;
    void calc_u(void);
    void modify_f(void);
    void calc_v(void);
    void calc_center(void);
    void export_conf(void);

    void print_positon(void);
};

#endif