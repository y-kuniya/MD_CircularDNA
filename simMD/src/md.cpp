#include "../include/md.h"
#include<random>
#include<iostream>
#include<fstream>

MD::MD(void){
    vars        = new Variables();
    vars_pred   = new Variables();
    obs         = new Observer();  
    force       = new Force();
}

MD::~MD(void){
    delete vars;
    delete vars_pred;
    delete obs;
    delete force;
}


// -----------------------初期配置の設定-----------------------------
// 粒子数を設定する関数
// makeconf()の前に必ず呼び出さなければいけない
void 
MD::set_PN(int PN_){
    PN  = PN_;
    SIZE= 3*PN;

    // beadの数が決まったのでここでD,D_sqrt,Fを確保して初期化する
    D_pred.resize(SIZE, vector<double>(SIZE,0.0));
    D_corr.resize(SIZE, vector<double>(SIZE,0.0));
    D_sqrt.resize(SIZE, vector<double>(SIZE,0.0));
    S.resize(SIZE, 0.0);
    vars->init(PN);
    vars_pred->init(PN);
    force->init(PN);
}

// 初期配置の設定 -- vars のr,b,uを設定する --
void
MD::makeconf(void){
    double r_sub[3];
    const double Beta = 2.0*M_PI/(double)PN;
    const double cosBeta = cos(Beta);
    const double sinBeta = sin(Beta);

    double *r = vars->r.data();
    double *f = vars->f.data();
    double *v = vars->v.data();
    double *u = vars->u.data();

    // 0番目の粒子をセット
    for(int i=0;i<3;i++){
        r[i] = 0.0;
    }
    
    f[0] = 1;
    f[1] = 0;
    f[2] = 0;

    v[0] = 0;
    v[1] = 1;
    v[2] = 0;

    u[0] = 0;
    u[1] = 0;
    u[2] = 1;

    for(int p=1;p<PN;p++){
        // r
        for(int i=0;i<3;i++){
            r[3*p+i] = r[3*(p-1)+i] + u[3*(p-1)+i]; 
        }
        // f
        for(int i=0;i<3;i++){
            f[3*p+i] = cosBeta*f[3*(p-1)+i] - sinBeta*u[3*(p-1)+i];
        }
        // u
        for(int i=0;i<3;i++){
            u[3*p+i] = sinBeta*f[3*(p-1)+i] + cosBeta*u[3*(p-1)+i];
        }
    }

    // b.inv_bの計算
    vars->calc_u();
    
}


// -----------------------予測子の計算-----------------------------
void 
MD::update_pred(vector<double> &R){
    double *r = vars->r.data();
    double *F = force->F.data();
    double *r_pred = vars_pred->r.data();

    // 1. 流体力学的テンソルの計算
    make_D(D_pred,vars->r,PN);
    // 2. 力の計算
    force->calc_force(vars);
    // 3. Sの計算,Sによるvars_predの更新
    for(int i=0;i<SIZE;i++){
        S[i] = 0.0;
        for(int j=0;j<SIZE;j++){
            S[i]+= D_pred[i][j]*F[j];
        }
        r_pred[i] = r[i]+S[i]*dt;
    }
    // 4. 揺動力による位置の更新
    // コレスキー分解
    Cholesky_decomposition(D_sqrt,D_pred,SIZE);
    // 確率変数を変換して、位置を更新
    for(int i=0;i<SIZE;i++){
        for(int j=0;j<SIZE;j++){
            r_pred[i] += D_sqrt[i][j]*R[j];
        }
    }
    // 5. r_predの更新が終わったので、フレームの計算をする
    vars_pred->calc_u();
}

// -----------------------修正子の計算-----------------------------
void 
MD::update_corr(vector<double> &R){
    double *r = vars->r.data();
    double *F = force->F.data();

    // 1. 流体力学的テンソルの計算
    make_D(D_corr,vars_pred->r,PN);
    // 2.力の計算
    force->calc_force(vars_pred);
    // 3. 外力による位置の更新
    for(int i=0;i<SIZE;i++){
        for(int j=0;j<SIZE;j++){
            S[i]+= D_corr[i][j]*F[j];
        }
        r[i]+=0.5*S[i]*dt;
    }
    // 4. 揺動力による位置の更新
    // D_meanの計算
    for(int i=0;i<SIZE;i++){
        for(int j=0;j<SIZE;j++){
            D_pred[i][j] += D_corr[i][j];
            D_pred[i][j] *= 0.5;
        }
    }
    // コレスキー分解
    Cholesky_decomposition(D_sqrt,D_pred,SIZE);
    // 確率変数を変換して、位置を更新
    for(int i=0;i<SIZE;i++){
        for(int j=0;j<SIZE;j++){
            r[i] += D_sqrt[i][j]*R[j];
        }
    }
    // 5. 位置の計算が終わったのでフレームの計算
    vars->calc_u();
}


//-----------------------1つのtrajectoryを作り、confファイルを吐く-----------------------------
void 
MD::make_trajectory(void){
    const int STEP      = 10000;
    const int DISPOSE   = 0;
    const int INTERVAL  = 10;

    random_device rnd;
    mt19937 mt(rnd());
    normal_distribution<> dist(0.0,sqrt(2.0*dt));
    
    // 1. 初期配置の設定
    set_PN(20);
    makeconf(); // set_PN()でSIZEが確定してから呼ぶ
    vector<double> R(SIZE,0.0); // set_PN()でSIZEが確定してから呼ぶ

    // 2. 前処理　--平衡状態の生成--
    for(int step=0;step<DISPOSE;step++){
        // predetor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_pred(R);
        // corretor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_corr(R);
    }

    // 3. 時間発展
    for(int step=0;step<STEP;step++){
        // conf ファイルを書き出す。
        if (step%INTERVAL ==0){
             vars->export_conf();
            //  cout<<step<<endl;
            //  vars->print_positon();
        }
        if(step%1000==0){
            cout<<step<<endl;
        }

        // predetor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_pred(R);
        // corretor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_corr(R);
    }
}


// ----------------------squared distance of center of mass の計算------------------------------
void 
MD::calc_sd(map<double,double> &m_sd,
            const int trial,
            const int STEP,
            const int DISPOSE,
            const int INTERVAL){
    random_device rnd;
    mt19937 mt(rnd());
    normal_distribution<> dist(0.0,sqrt(2.0*dt));

    // 1. 初期配置の設定
    set_PN(10);
    makeconf();
    vector<double> R(SIZE,0.0);

    // 2. 前処理　--平衡状態の生成--
    for(int step=0;step<DISPOSE;step++){
        // predetor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_pred(R);
        // corretor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_corr(R);
    }

    // 3. 起点となる重心位置
    obs->set_c0(vars);

    // 3. 時間発展
    for(int step=0;step<STEP;step++){
        // predetor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_pred(R);
        // corretor
        for(int i=0;i<SIZE;i++){
            R[i] = dist(mt);
        }
        update_corr(R);
        // 平均2乗距離の計算
        if((step+1)%INTERVAL==0){
            double t = static_cast<double>(step+1)*dt;
            double sd=obs->squared_distance(vars);
            m_sd[t] = ((double)trial)*m_sd[t]+sd;
            m_sd[t]/= (double)(trial+1);
        }
    }
}   

void run_m_sd(void){
    const int STEP      = 5000;
    const int DISPOSE   = 1000;
    const int INTERVAL  = 100;
    const int TRIAL     = 200;

    map<double,double> m_sd;

     for(int trial =0;trial<TRIAL;trial++){
        if(trial%10==0){
            cout<<trial<<endl;
        } 
        MD *md;
        md = new MD();
        md->calc_sd(m_sd,trial,STEP,DISPOSE,INTERVAL);
    }

    // 平均２乗距離をファイルへ書き出す
    string filename("../data/m_sd.csv");
    fstream ofs;
    ofs.open(filename, std::ios_base::out);
    for(auto &m:m_sd){
        ofs<<m.first<<","<<m.second<<endl;
    }
    ofs.close();
}