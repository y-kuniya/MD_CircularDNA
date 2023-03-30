#include "../include/force.h"
#include<cmath>

void
Force::init(const int PN_){
    PN  = PN_;
    SIZE= 3*PN;
    F.resize(SIZE,0.0);
    cosBeta.resize(PN-1,0.0);
    Beta_over_sinBeta.resize(PN-1,0.0);
}

void 
Force::clear(void){
    for(auto &f:F){
        f = 0.0;
    }
}

// ------------------------------------力の計算-------------------------------------------
void 
Force::calc_force(Variables *vars){
    // 力を0.0にセットする。
    clear();
    // 排除体積力の計算
    calc_excluedVolume(vars);
    // 伸び弾性力の計算
    calc_stretching(vars);
    // 曲げ角の計算
    calc_Beta(vars);
    // 曲げ弾性力の計算
    calc_bending(vars);
}

// #include<iostream>
// ----------------------------------excluded volume-------------------------------------
void 
Force::calc_excluedVolume(Variables *vars){
    double eps_over_cutoff = epsilon/cutoff;

    double dr[3];
    double r_pq,ratio,ratio_6,ratio_7;
    double Lpq,val;

    double *r = vars->r.data();

    for(int p=0;p<PN-1;p++){
        for(int q=p+1;q<PN;q++){
            // dr,r_pqの計算
            r_pq = 0.0;
            for(int i=0;i<3;i++){
                dr[i] = r[3*p+i] - r[3*q+i];
                r_pq += dr[i]*dr[i]; 
            }
            r_pq = sqrt(r_pq);

            // 排除体積が働く範囲か判定
            if (r_pq >= cutoff) continue;
            // 排除体積力の計算
            ratio = cutoff / r_pq ;
            ratio_6 = pow(ratio,6);
            ratio_7 = ratio_6*ratio;

            Lpq = 6*eps_over_cutoff*ratio_7*(2.0*ratio_6-1.0);
            Lpq /= r_pq;
            for(int i=0;i<3;i++){
               val = Lpq*dr[i];
               F[3*p+i] += val;
               F[3*q+i] -= val;
            }
        }
    }
}

// ----------------------------------stretching------------------------------------------
void 
Force::calc_stretching(Variables *vars){
    const double one_over_delta_delta = 1.0/(delta*delta);
    double coeff,val;

    double *u = vars->u.data();
    double *b = vars->b.data();

    for(int p=0;p<PN-1;p++){
        coeff = (b[p]-l_eq)*one_over_delta_delta;
        for(int i=0;i<3;i++){
            val = coeff*u[3*p+i];
            F[3*p+i]    += val ;
            F[3*(p+1)+i]-= val ;
        }
    }
}

// ----------------------------------Bending------------------------------------------
void
Force::calc_Beta(Variables *vars){
    double *u = vars->u.data();

    // cos(Beta)
    for(int p=0;p<PN-1;p++){
        cosBeta[p] = 0.0;
        for(int i=0;i<3;i++){
            cosBeta[p] += u[3*(p+1)+i]*u[3*p+i];
        }
    }

    // Beta / sin(Beta)
    double angle;
    for(int p=0;p<PN-1;p++){
        angle = acos(cosBeta[p]);
        if (angle <1e-4){
            Beta_over_sinBeta[p] = 1.0;
        }else{
            Beta_over_sinBeta[p] = angle/sqrt(1.0 - cosBeta[p]*cosBeta[p]);
        }
    }
}


void 
Force::calc_bending(Variables *vars){
    double one_over_psi_psi = 1.0/(psi*psi);
    double val; 

    double *u = vars->u.data();
    double *inv_b = vars->inv_b.data();

     // Aについての更新
    for(int p=0;p<PN-2;p++){
        for(int i=0;i<3;i++){
            val = one_over_psi_psi*Beta_over_sinBeta[p]*inv_b[p]*( u[3*(p+1)+i]-u[3*p+i]*cosBeta[p] );
            F[3*p+i]    += -val;
            F[3*(p+1)+i]+=  val;
        }
    }

    // Bについての更新
    for(int p=1;p<PN-1;p++){
        for(int i=0;i<3;i++){
            val = one_over_psi_psi*Beta_over_sinBeta[p-1]*inv_b[p]*(u[3*(p-1)+i]-u[3*p+i]*cosBeta[p-1]);
            F[3*p+i]    += -val;
            F[3*(p+1)+i]+=  val;
        }
    }
}