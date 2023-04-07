#include "../include/tools.h"
#include "../include/systemparam.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void print_matrix(const SquareVec &a,const int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<< scientific << setprecision(2) << uppercase <<a[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl;
}

// 行列Bを行列Aに掛ける
void product_matrix(SquareVec &A, const SquareVec &B){
    SquareVec C(3,vector<double>(3,0.0));

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            A[i][j] = C[i][j];
        }
    }
}

// 一番最初に行列Rが0.0に初期化されている必要がある
void set_R_3(SquareVec &Rotation,const double angle){
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);

    Rotation[0][0] = cosAngle;
    Rotation[0][1] = -sinAngle;
    Rotation[1][0] = sinAngle;
    Rotation[1][1] = cosAngle;
    Rotation[2][2] = 1.0;
}

void set_R_2(SquareVec &Rotation,const double angle){
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);

    Rotation[0][0] = cosAngle;
    Rotation[0][2] = sinAngle;
    Rotation[1][1] = 1.0;
    Rotation[2][0] = -sinAngle;
    Rotation[2][2] = cosAngle;
}

void Cholesky_decomposition(SquareVec &L,const SquareVec &a,const int n){
    double s;
    for(int i=0;i<n;i++){
        for(int j=0;j<i;j++){
            s = a[i][j];
            for(int k=0;k<j;k++) s -= L[i][k] * L[j][k];
            L[i][j] = s / L[j][j];
        }
        /* L[i][i] の計算 */
        s = a[i][i];
        for(int k=0;k<i;k++) s -= L[i][k]*L[i][k];

        if (s < 0){
            cout<<"s < 0"<<endl;
            exit(0);
        }

        L[i][i] = sqrt(s);
    }
}

// 流体力学的テンソルの計算
// ずっと0のところにはアクセスしないので配列を確保したとき0.0に初期化してあることが前提
void make_D(SquareVec &D,vector<double> &r,const int PN){
    double sigma_sigma_2 = 2.0*sigma*sigma;
    double three_over_sigma_32 = 3.0/(32.0*sigma);
    double nine_over_sigma_32  = 9.0/(32.0*sigma);

    double delta_r[3];
    double coeff,r_pq,r_pq_2;
    double term_1,term_2;

    for(int p = 0;p<PN;p++){
        // 非対角成分
        for(int q=0;q<p;q++){
            // ---前処理---
            r_pq_2 = 0.0;
            for(int i=0;i<3;i++){
                delta_r[i] = r[3*p+i] - r[3*q+i];
                r_pq_2 += delta_r[i]*delta_r[i];
            }
            r_pq = sqrt(r_pq_2);
            
            // 粒子が重ならない場合
            if(r_pq > l_eq){
                coeff= (3.0*sigma)/(4.0*r_pq);
                for(int i=0;i<3;i++){
                    for(int j=0;j<3;j++){
                        term_1 = delta_r[i]*delta_r[j]/r_pq_2;
                        term_2 = -3.0*term_1;
                        
                        if(i==j){
                            term_1 += 1.0;
                            term_2 += 1.0;
                        }
                        term_2 *= sigma_sigma_2/(3.0*r_pq_2);
                        
                        D[3*p+i][3*q+j] = coeff*(term_1 + term_2);
                        // 対称部分にも代入
                        D[3*q+j][3*p+i] = D[3*p+i][3*q+j];
                    }
                }
            }
            // 粒子が重なる場合
            else{
                for(int i=0;i<3;i++){
                    for(int j=0;j<3;j++){
                        term_1 = delta_r[i]*delta_r[j]/r_pq_2;
                        term_1*= three_over_sigma_32;
                        if(i==j){
                            term_1 += 1.0 - r_pq*nine_over_sigma_32;
                        }

                        D[3*p+i][3*q+j] = term_1;
                        // 対称部分にも代入
                        D[3*q+j][3*p+i] = D[3*p+i][3*q+j];
                    }
                }
            }
        }

        // 対角成分(p==q)
        // ブロック行列の対角成分 [3*p+i][3*p+i]にだけアクセス
        for(int i=0;i<3;i++){
            D[3*p+i][3*p+i] = 1.0;
        }
    }
}