#include "../include/observer.h"
#include <cmath>

Observer::Observer(void){
    string filename("../data/topologicalPropaties.csv");
    ofs.open(filename, std::ios_base::out);
}

Observer::~Observer(void){
    ofs.close();
}

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

double 
Observer::writhing_number(Variables *vars){
    double *r = vars->r.data();
    double *u = vars->u.data();
    double *b = vars->b.data();
    const int PN = vars->num;

    double r_pq[3];
    double bpxbq[3];
    double norm_r_pq;
    double inv_norm_r_pq;
    double inv_3;
    

    double Wr = 0.0;
    double coeff;
    // for(int p=0;p<PN-1;p++){
    //     for(int q=p+1;q<PN;q++){
    
    for(int p=0;p<PN;p++){
        for(int q=0;q<PN;q++){
            if(p==q) continue;

            norm_r_pq = 0.0;
            for(int i=0;i<3;i++){
                r_pq[i] = r[3*p+i] - r[3*q+i];
                norm_r_pq += r_pq[i]*r_pq[i];
            }
            norm_r_pq = sqrt(norm_r_pq);
            inv_norm_r_pq = 1.0/norm_r_pq;
            inv_3 = inv_norm_r_pq*inv_norm_r_pq*inv_norm_r_pq; 

            bpxbq[0] = u[3*p+1]*u[3*q+2] - u[3*p+2]*u[3*q+1];
            bpxbq[1] = u[3*p+2]*u[3*q] - u[3*p]*u[3*q+2];
            bpxbq[2] = u[3*p]*u[3*q+1] - u[3*p+1]*u[3*q];

            coeff = b[p]*b[q]; 
            for(int i=0;i<3;i++){
                Wr += coeff*bpxbq[i]*r_pq[i]*inv_3;
            }
        }
    }

    // Wr /= 2.0*M_PI;
    Wr /= 4.0*M_PI;
    return Wr;
}

double 
Observer::total_twist(Force *force){
    double Tw = 0.0;

    for(auto& ag:force->AG){
        Tw += ag;
    }
    Tw/= 2.0*M_PI;
    
    return Tw;
}

void
Observer::export_Wr_Tw(int step,Variables *vars,Force *force){
    double Wr = writhing_number(vars);
    double Tw = total_twist(force);
    
    // double DLk;
    // if(Lk==0)   DLk = 100.0*(Wr+Tw-Lk);
    // else        DLk= 100.0*(Wr+Tw-Lk)/Lk;
    
    ofs<<(step+1)<<","<<Wr<<","<<Tw<<","<<Wr+Tw<<endl;
}