#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <vector>
using SquareVec = std::vector<std::vector<double>>;

void print_matrix(const SquareVec &a,const int n);
void Cholesky_decomposition(SquareVec &L,const SquareVec &a,const int n);
void make_D(SquareVec &D,std::vector<double> &r,const int PN);

#endif 