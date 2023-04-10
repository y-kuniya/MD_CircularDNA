#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <vector>
using SquareVec = std::vector<std::vector<double>>;

void print_matrix(const SquareVec &a,const int n);
void product_matrix(SquareVec &A,const SquareVec &B);
void set_R_3(SquareVec &Rotation,const double angle);
void set_R_2(SquareVec &Rotation,const double angle);
void Cholesky_decomposition(SquareVec &L,const SquareVec &a,const int n);
void make_D(SquareVec &D,std::vector<double> &r,const int PN);

#endif 