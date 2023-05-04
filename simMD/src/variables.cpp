#include "../include/variables.h"
#include <fstream>
#include <cmath> // sqrt at calc_u()
 
Variables::Variables(void){
    num  = 0 ;
}

void 
Variables::init(int num_){
    num = num_;
    const int SIZE = 3*num;

    r.resize(SIZE,0.0);
    f.resize(SIZE,0.0);
    v.resize(SIZE,0.0);
    u.resize(SIZE,0.0);
    b.resize(num,0.0);
    inv_b.resize(num,0.0);
}

void 
Variables::calc_u(void){
    const int PN = num;
    for(int p=0;p<PN;p++){
        double sd = 0.0; // squared distance
        for(int i=0;i<3;i++){   
            // 端
            if (p==PN-1){
                u[3*p+i] = r[i] - r[3*p+i];
            // バルク
            }else{
                u[3*p+i] = r[3*(p+1)+i] - r[3*p+i];
            }
            sd += u[3*p+i]*u[3*p+i];
        }
        // bを格納
        b[p] = sqrt(sd);
        // 規格化
        inv_b[p] = 1.0/b[p];
        for(int i=0;i<3;i++){
            u[3*p+i] *= inv_b[p];
        }
    }
}

void 
Variables::modify_f(void){
    const int PN = num;
    double f_dot_u,norm_f,inv_norm_f;

    for(int p=0;p<PN;p++){
        // fからuの方向を引く
        f_dot_u = 0.0;
        for(int i=0;i<3;i++){
            f_dot_u += f[3*p+i]*u[3*p+i];
        }
        for(int i=0;i<3;i++){
            f[3*p+i] -= f_dot_u*u[3*p+i];
        }
        // fを規格化する
        norm_f = 0.0;
        for(int i=0;i<3;i++){
            norm_f += f[3*p+i]*f[3*p+i];
        }
        norm_f = sqrt(norm_f);
        inv_norm_f = 1.0/norm_f;
        for(int i=0;i<3;i++){
            f[3*p+i] *= inv_norm_f;
        }
    }
}

void 
Variables::calc_v(void){
    const int PN = num;
    for(int p=0;p<PN;p++){
        v[3*p]     = u[3*p+1]*f[3*p+2] - u[3*p+2]*f[3*p+1];
        v[3*p+1]   = u[3*p+2]*f[3*p] - u[3*p]*f[3*p+2];
        v[3*p+2]   = u[3*p]*f[3*p+1] - u[3*p+1]*f[3*p]; 
    }   
}

void 
Variables::calc_center(void){
    for(int i=0;i<3;i++){
        c[i] = 0.0;
    }
    for(int p=0;p<num;p++){
        for(int i=0;i<3;i++){
            c[i] += r[3*p+i];
        }
    }

    for(int i=0;i<3;i++){
        c[i] /= static_cast<double>(num);
    }
} 

void 
Variables::export_conf(void){
    static int count = 0;
    char filename[256];
    sprintf(filename,"../confData/conf%03d.bin",count);
    count++;

    ofstream ofs(filename, std::ios::out | std::ios::binary);
    ofs.write((char*)&r[0], r.size()*sizeof(r[0]));
    ofs.close();
}

void 
Variables::export_conf_csv(void){

    char filename[256];
    sprintf(filename,"../data/conf.csv");
    ofstream ofs(filename, std::ios::out);

    calc_center();
    for(int p = 0;p<num;p++){
        ofs<<r[3*p] - c[0]<<","<<r[3*p+1] - c[1]<<","<<r[3*p+2] - c[2]<<endl;
    }
}

#include <iostream>
void 
Variables::print_positon(void){
    for(int p=0;p<num;p++){
        cout<<p<<" | ";
        for(int i=0;i<3;i++){
            cout<<r[3*p+i]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}