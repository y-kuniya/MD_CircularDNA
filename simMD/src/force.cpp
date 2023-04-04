#include "../include/force.h"
#include<cmath>

void
Force::init(const int PN_){
    PN  = PN_;
    SIZE= 3*PN;

    F.resize(SIZE,0.0);
    // 曲げ角
    cosBeta.resize(PN,0.0);
    Beta_over_sinBeta.resize(PN,0.0);
    // ねじれ角
    cosAG.resize(PN,0.0);
    sinAG.resize(PN,0.0);
    AG.resize(PN,0.0);
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

    // 周期境界の所の計算
    coeff = (b[PN-1]-l_eq)*one_over_delta_delta;
    for(int i=0;i<3;i++){
        val = coeff*u[3*(PN-1)+i];
        F[3*(PN-1)+i]   += val;
        F[i]            -= val;
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

    // 周期境界条件
    cosBeta[PN-1]=0.0;
    for(int i=0;i<3;i++){
        cosBeta[PN-1] += u[i]*u[3*(PN-1)+i];
    } 

    // Beta / sin(Beta)
    double angle;
    for(int p=0;p<PN;p++){
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
    for(int p=0;p<PN-1;p++){
        for(int i=0;i<3;i++){
            val = one_over_psi_psi*Beta_over_sinBeta[p]*inv_b[p]*( u[3*(p+1)+i]-u[3*p+i]*cosBeta[p] );
            F[3*p+i]    += -val;
            F[3*(p+1)+i]+=  val;
        }
    }

    for(int i=0;i<3;i++){
            val = one_over_psi_psi*Beta_over_sinBeta[PN-1]*inv_b[PN-1]*( u[i]-u[3*(PN-1)+i]*cosBeta[PN-1] );
            F[3*(PN-1)+i]   += -val;
            F[i]            +=  val;
    }

    
    // Bについての更新
    for(int i=0;i<3;i++){
            val = one_over_psi_psi*Beta_over_sinBeta[PN-1]*inv_b[0]*(u[3*(PN-1)+i]-u[i]*cosBeta[PN-1]);
            F[i]    += -val;
            F[3+i]  +=  val;
    }

    for(int p=1;p<PN;p++){
        for(int i=0;i<3;i++){
            val = one_over_psi_psi*Beta_over_sinBeta[p-1]*inv_b[p]*(u[3*(p-1)+i]-u[3*p+i]*cosBeta[p-1]);
            F[3*p+i]    += -val;
            F[3*(p+1)+i]+=  val;
        }
    }
}


// ----------------------------------Torsional------------------------------------------
// 必ずBetaの計算の後に呼び出すこと
void 
Force::calc_AlphaPlusGamma(Variables *vars){
    double *f = vars->f.data();
    double *v = vars->v.data();
    double *u = vars->u.data();

    // cos(Alpha + Gamma) 
    for(int p=0;p<PN-1;p++){
        cosAG[p] = 0.0;
        for(int i=0;i<3;i++){
            cosAG[p] += f[3*p+i]*f[3*(p+1)+i] + v[3*p+i]*v[3*(p+1)+i]; 
        }
        cosAG[p] /= (1.0 + cosBeta[p]);
    }
    
    cosAG[PN-1] = 0.0;
    for(int i=0;i<3;i++){
        cosAG[PN-1] += f[3*(PN-1)+i]*f[i] + v[3*(PN-1)+i]*v[i];
    }
    cosAG[PN-1] /= (1.0+cosBeta[PN-1]);

    // sin(Alpha + Gamma)
    for(int p=0;p<PN-1;p++){
        sinAG[p] = 0.0;
        for(int i=0;i<3;i++){
            sinAG[p] += v[3*p+i]*f[3*(p+1)+i] - f[3*p+i]*v[3*(p+1)+i];
        }
        sinAG[p] /= (1.0+cosBeta[p]);
    }

    sinAG[PN-1] = 0.0;
    for(int i=0;i<3;i++){
        sinAG[PN-1] += v[3*(PN-1)+i]*f[i] - f[3*(PN-1)+i]*v[i];
    }
    sinAG[PN-1] /= (1.0+cosBeta[PN-1]);

    // Alpha + Gamma 
    for(int p=0;p<PN;p++){
        AG[p] = acos(cosAG[p]);
        if (sinAG[p] < 0.0) AG[p] *= -1.0;
    }
}


void
Force::calc_torsional(Variables *vars){

    double *f = vars->f.data();
    double *v = vars->v.data();
    double *u = vars->u.data();
    double *inv_b = vars->inv_b.data();

}