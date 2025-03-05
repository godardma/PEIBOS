// Author: Maël GODARD

#include "peibos/PEIBOS.h"


int main()
{
  VectorVar y(3);
  AnalyticFunction f({y},{sqr(y[0])-sqr(y[1])+y[0],2*y[0]*y[1]+y[1],y[2]});
  AnalyticFunction f_id({y},{y[0],y[1],y[2]});
  
  VectorVar X(2);
  AnalyticFunction psi0 ({X},{1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

  vector<initializer_list<int>> generators ({{1,2,3},
                                            {-2,1,3},
                                            {3,2,-1}});

  double epsilon = 0.1;
  
  PEIBOS(f_id, psi0, generators, epsilon, "Atlas");
  PEIBOS(f, psi0, generators, epsilon, "Conform");
}

  
