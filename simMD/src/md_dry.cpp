#include "../include/md_dry.h"
#include<random>
#include<iostream>
#include<fstream>

MDdry::MDdry(void){
    vars        = new Variables();
    vars_pred   = new Variables();
    obs         = new Observer();  
    force       = new Force();
}

MDdry::~MDdry(void){
    delete vars;
    delete vars_pred;
    delete obs;
    delete force;
}

// -----------------------初期配置の設定-----------------------------
// 粒子数を設定する関数
// makeconf()の前に必ず呼び出さなければいけない
void 
MDdry::set_PN(int PN_){
    PN  = PN_;
    SIZE= 3*PN;

    // beadの数が決まったのでここでD,D_sqrt,Fを確保して初期化する
    S.resize(SIZE, 0.0);
    DT.resize(SIZE,0.0);
    u_sub.resize(SIZE,0.0);
    vars->init(PN);
    vars_pred->init(PN);
    force->init(PN);
}

// 初期値の設定  varsのr,f,v,uを設定する
void 
MDdry::makeconf(void){
    const double delta_Lk = -12.0;
    // const double delta_Lk = 0.0;
    const double delta_twist = (2.0*M_PI)*delta_Lk/(double)PN; 

    double *r = vars->r.data();
    double *f = vars->f.data();
    double *v = vars->v.data();
    double *u = vars->u.data();
   
    SquareVec Frame(3,vector<double>(3,0.0));
    SquareVec R_alpha(3,vector<double>(3,0.0));
    SquareVec R_beta(3,vector<double>(3,0.0));
    SquareVec R_gamma(3,vector<double>(3,0.0));

    double alpha,beta,gamma_angle;

    // ---1. フレームの設定---
    // y軸回転は一定なので最初に設定する
    beta = -2.0*M_PI/(double)PN;
    set_R_2(R_beta,beta);

    // 0番目のbeadのフレームを設定する
    for(int i=0;i<3;i++){
        Frame[i][i] = 1.0;
    }
    // フレームをvarsの変数に格納
    for(int i=0;i<3;i++){
        f[i] = Frame[i][0];
        v[i] = Frame[i][1];
        u[i] = Frame[i][2];
    }

    for(int p=1;p<PN;p++){
        // Frameの回転
        alpha = -(double)(p-1)*delta_twist;
        gamma_angle = (double)p*delta_twist;
        set_R_3(R_alpha,alpha);
        set_R_3(R_gamma,gamma_angle);

        product_matrix(Frame,R_alpha);
        product_matrix(Frame,R_beta);
        product_matrix(Frame,R_gamma);

        // Frameをvarsの変数に格納
        for(int i=0;i<3;i++){
        f[3*p+i] = Frame[i][0];
        v[3*p+i] = Frame[i][1];
        u[3*p+i] = Frame[i][2];
        }
    }

    // ---2. 位置を設定---
    for(int i=0;i<3;i++){
        r[i] = 0.0;
    }
    for(int p=1;p<PN;p++){
        for(int i=0;i<3;i++){
            r[3*p+i] = r[3*(p-1)+i] + u[3*(p-1)+i]; 
        }
    }

    // ---3.ボンド長の計算---
    vars->calc_u();
}

// -----------------------予測子の計算-----------------------------
void 
MDdry::update_pred(vector<double> &R){
    double *r = vars->r.data();
    double *F = force->F.data();
    double *r_pred = vars_pred->r.data();

    // 1. 力の計算
    force->calc_force(vars);
    // 2. Sの計算,Sによるvars_predの更新
    for(int i=0;i<SIZE;i++){
        S[i] = F[i];
        r_pred[i] = r[i]+S[i]*dt;
    }
    // 3. 揺動力による位置の更新
    for(int i=0;i<SIZE;i++){
         r_pred[i] += R[i];
    }
}

// -----------------------修正子の計算-----------------------------
void 
MDdry::update_corr(vector<double> &R){
    double *r = vars->r.data();
    double *F = force->F.data();

    // 1.力の計算
    force->calc_force(vars_pred);
    // 2. 外力による位置の更新
    for(int i=0;i<SIZE;i++){
        S[i] += F[i];
        r[i]+=0.5*S[i]*dt;
    }
    // 3. 揺動力による位置の更新
    for(int i=0;i<SIZE;i++){
        r[i] += R[i];
    }
}

// -----------------------フレームに関した予測子の計算-----------------------
void 
MDdry::update_frame_pred(vector<double> &R){
    double delta_phi,f_dot_u;
    double uxf[3];

    double *f = vars->f.data();
    double *u = vars->u.data();
    double *f_pred = vars_pred->f.data();
    double *u_pred = vars_pred->u.data();
    double *T = force->T.data();

    // 1. u_predの計算
    vars_pred->calc_u();
    // 2. ねじれ弾性によるトルクの計算
    force->calc_torque();
    // 3. fの更新
    for(int p=0;p<PN;p++){    
        // 2-1. delta_phiの計算
        DT[p] = Dr*T[p]*dt;
        delta_phi = DT[p] + R[p];
        // 2-2. f_predの更新
        f_dot_u = 0.0;
        for(int i=0;i<3;i++){
            f_dot_u += f[3*p+i]*u_pred[3*p+i];
        }
        uxf[0] = u[3*p+1]*f[3*p+2] - u[3*p+2]*f[3*p+1];
        uxf[1] = u[3*p+2]*f[3*p] - u[3*p]*f[3*p+2];
        uxf[2] = u[3*p]*f[3*p+1] - u[3*p+1]*f[3*p];
        
        for(int i=0;i<3;i++){
            f_pred[3*p+i] = f[3*p+i] + delta_phi*uxf[i] - u[3*p+i]*f_dot_u; 
        }
    }
    // 3. f_predの修正
    vars_pred->modify_f();
    // 4. v_predの計算
    vars_pred->calc_v();
}

// -----------------------フレームに関した修正子の計算-----------------------
void 
MDdry::update_frame_corr(vector<double> &R){
    double delta_phi,f_dot_u;
    double uxf[3];

    double *f = vars->f.data();
    double *u = vars->u.data();
    double *T = force->T.data();

    // 1. uの更新
    for(int pi=0;pi<SIZE;pi++){
        u_sub[pi] = u[pi];
    }
    vars->calc_u();
    // 2. ねじれ弾性によるトルクの計算
    force->calc_torque();
    // 3.fの更新
    for(int p=0;p<PN;p++){
        // 2-1. delta_phiの計算
        DT[p] += Dr*T[p]*dt;
        DT[p] *= 0.5;
        delta_phi = DT[p] + R[p];
        // 2-2. fの更新
        f_dot_u = 0.0;
        for(int i=0;i<3;i++){
            f_dot_u += f[3*p+i]*u[3*p+i];
        }
        uxf[0] = u_sub[3*p+1]*f[3*p+2] - u_sub[3*p+2]*f[3*p+1];
        uxf[1] = u_sub[3*p+2]*f[3*p] - u_sub[3*p]*f[3*p+2];
        uxf[2] = u_sub[3*p]*f[3*p+1] - u_sub[3*p+1]*f[3*p];
        
        for(int i=0;i<3;i++){
            f[3*p+i] = f[3*p+i] + delta_phi*uxf[i] - u_sub[3*p+i]*f_dot_u; 
        }
    }
    // 4. fの修正
    vars->modify_f();
    // 5. vの計算
    vars->calc_v();
}

//-----------------------1つのtrajectoryを作り、confファイルを吐く-----------------------------
void 
MDdry::make_trajectory(void){
    const int STEP      = 400000;
    const int DISPOSE   = 0;
    const int INTERVAL  = 200;

    random_device rnd;
    mt19937 mt(rnd());
    normal_distribution<> dist(0.0,sqrt(2.0*dt));
    normal_distribution<> dist_phi(0.0,sqrt(2.0*Dr*dt));

    // 1. 初期配置の設定
    set_PN(60);
    makeconf(); // set_PN()でSIZEが確定してから呼ぶ
    vector<double> R(SIZE,0.0); // set_PN()でSIZEが確定してから呼ぶ
    vector<double> R_phi(PN,0.0);

    // 3. 時間発展
    for(int step=0;step<STEP;step++){
        // conf ファイルを書き出す。
        if (step%INTERVAL ==0){
             vars->export_conf();
            //  vars->print_positon();
            //  cout<<step<<endl;
        }
        if(step%5000==0){
            cout<<step<<" : ";
            double Wr = obs->writhing_number(vars);
            double Tw = force->calc_total_twist();
            cout<<-Wr<<" ";
            // cout<< Tw <<endl;
            cout<<100*(Wr+Tw+12.0)/12.0 <<endl;
        }


        // ---predetor--- 
        // 位置
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_pred(R);
        // frame
        for(int i=0;i<PN;i++){
            R_phi[i] = dist_phi(mt);
        }
        update_frame_pred(R_phi);

        // ---corretor---
        // 位置
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_corr(R);
        // frame
        for(int i=0;i<PN;i++){
            R_phi[i] = dist_phi(mt);
        }
        update_frame_corr(R_phi);
    }

}

